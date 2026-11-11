! Floating-default maximum-entropy Doppler tomography: numerical core.
!
! Based on H. C. Spruit's 1999 code (dopmap, version 2.0). Revised forward
! model, likelihoods, optimizer, regularization schedule and storage.
! All arrays are allocated at run time, so one executable serves every
! problem size. Data index k = ip + (iv-1)*nph (phase varies fastest).
module dop_science
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_fortran_env, only: int32, real64
  implicit none

  integer, parameter :: dp = real64
  real(dp), parameter :: pi2 = 8.0_dp*atan(1.0_dp)
  integer, parameter :: max_half_width = 100

  ! Problem dimensions.
  integer :: nph = 0, nvp = 0, nv = 0, nil = 0, nol = 0
  ! Statistical and physical model.
  integer :: stat_iw = 0, phase_subsamples = 0, statistic_rank = 0
  integer :: smoothing_width = 7, log_unit = 0
  ! Exposure-quadrature samples per (smaller of map pixel, channel) of
  ! radial-velocity travel; 3 keeps errors below ~0.05% of a pixel's peak response.
  real(dp) :: quadrature_density = 3
  real(dp) :: flux_scale = 1, instrument_sigma = 0, baseline_fraction = 0.1_dp
  real(dp) :: baseline_level = 0, likelihood_scale = 1, information = 1
  real(dp) :: solver_objective = 0, solver_stationarity = 0
  real(dp) :: wid = 1, af = 0, va = 0, dvmap = 0
  logical :: baseline_raised = .false.
  real(dp) :: rms_floor_estimate = -1
  real(dp), allocatable :: precision(:,:,:), stat_weights(:), normalization(:), dpha(:)
  integer, allocatable :: covariance_rank(:), indx(:), indy(:)
  logical, allocatable :: phase_valid(:)
  ! Pixel-major sparse projection: pixel j contributes proj_value(k) to datum
  ! proj_index(k) for k = proj_ptr(j), ..., proj_ptr(j+1)-1.
  integer, allocatable :: proj_ptr(:), proj_index(:)
  real(dp), allocatable :: proj_value(:)
  ! Last four gradient maps, reported as diagnostic "change maps".
  real, allocatable :: dpx(:,:)
  ! Cached default-map smoothing kernel.
  integer :: kernel_ns = -1, kernel_half = 0
  real(dp) :: kernel(-max_half_width:max_half_width) = 0

contains

  subroutine fail(message)
    character(len=*), intent(in) :: message
    write(*,'(a)') trim(message)
    if (log_unit /= 0) write(log_unit,'(a)') trim(message)
    stop 1
  end subroutine fail

  subroutine apply_precision(residual, weighted, quadratic)
    ! weighted = W r per spectrum (W diagonal or full precision); masked
    ! spectra contribute nothing. quadratic = r^T W r.
    real(dp), intent(in) :: residual(:)
    real(dp), intent(out) :: weighted(:), quadratic
    real(dp) :: rv(nvp), zv(nvp)
    integer :: ip, iv
    quadratic = 0
    weighted(1:nil) = 0
    do ip = 1, nph
      if (.not. phase_valid(ip)) cycle
      do iv = 1, nvp
        rv(iv) = residual(ip + (iv-1)*nph)
      end do
      if (stat_iw == 2) then
        zv = matmul(precision(:,:,ip), rv)
      else
        do iv = 1, nvp
          zv(iv) = stat_weights(ip + (iv-1)*nph)*rv(iv)
        end do
      end if
      do iv = 1, nvp
        weighted(ip + (iv-1)*nph) = zv(iv)
      end do
      quadratic = quadratic + dot_product(rv, zv)
    end do
  end subroutine apply_precision

  subroutine set_phase_mask(pha, pb0, pb1)
    ! Open, periodically repeated exclusion interval (pb0, pb1) in radians.
    real(dp), intent(in) :: pha(:), pb0, pb1
    real(dp) :: width, tolerance, distance
    integer :: i
    if (allocated(phase_valid)) deallocate(phase_valid)
    allocate(phase_valid(nph))
    width = (pb1 - pb0)/pi2
    tolerance = 16*epsilon(1.0_dp)
    do i = 1, nph
      distance = modulo((pha(i) - pb0)/pi2, 1.0_dp)
      phase_valid(i) = .not. (distance > tolerance .and. distance < width - tolerance)
    end do
  end subroutine set_phase_mask

  subroutine statistical_input(dat, eb, iw, norm)
    real(dp), intent(inout) :: dat(:), eb(:)
    integer, intent(in) :: iw, norm
    integer :: ip, iv, k, u, ios
    integer(int32) :: n1, n2, rank
    real(dp) :: avg, su, vmax
    if (allocated(normalization)) deallocate(normalization)
    if (allocated(covariance_rank)) deallocate(covariance_rank)
    allocate(normalization(nph), covariance_rank(nph))
    normalization = 1
    covariance_rank = nvp
    if (count(phase_valid) == 0) call fail('No unmasked orbital phases remain.')
    if (iw < 0 .or. iw > 2) call fail('Invalid weight mode iw.')
    do k = 1, nil
      if (.not. ieee_is_finite(dat(k))) call fail('Non-finite input data.')
      if (iw >= 1) then
        if (.not. ieee_is_finite(eb(k)) .or. eb(k) <= 0) call fail('Invalid error bars.')
      else
        eb(k) = 1
      end if
    end do
    if (iw == 2) then
      if (allocated(precision)) deallocate(precision)
      allocate(precision(nvp, nvp, nph))
      open(newunit=u, file='dop.precision', access='stream', form='unformatted', &
           status='old', action='read', convert='little_endian', iostat=ios)
      if (ios /= 0) call fail('Cannot open dop.precision.')
      read(u, iostat=ios) n1, n2
      if (ios /= 0 .or. n1 /= nph .or. n2 /= nvp) call fail('dop.precision has the wrong shape.')
      do ip = 1, nph
        read(u, iostat=ios) rank
        if (ios /= 0 .or. rank < 1 .or. rank > nvp) call fail('Invalid covariance rank.')
        covariance_rank(ip) = rank
        read(u, iostat=ios) precision(:, :, ip)
        if (ios /= 0) call fail('Truncated dop.precision.')
        if (.not. all(ieee_is_finite(precision(:, :, ip)))) call fail('Non-finite precision.')
      end do
      close(u)
    end if
    ! Optional normalization to the mean continuum-subtracted line integral.
    if (norm == 1) then
      avg = 0
      do ip = 1, nph
        if (.not. phase_valid(ip)) cycle
        su = 0
        do iv = 1, nvp
          su = su + dat(ip + (iv-1)*nph)
        end do
        if (su <= 0) call fail('Phase normalization requires positive line flux.')
        normalization(ip) = su
        avg = avg + su
      end do
      avg = avg/count(phase_valid)
      do ip = 1, nph
        if (phase_valid(ip)) normalization(ip) = avg/normalization(ip)
      end do
    end if
    do ip = 1, nph
      do iv = 1, nvp
        k = ip + (iv-1)*nph
        dat(k) = dat(k)*normalization(ip)
        eb(k) = eb(k)*normalization(ip)
      end do
    end do
    ! An internal RMS flux scale makes all numerical choices unit invariant.
    vmax = 0
    do ip = 1, nph
      if (.not. phase_valid(ip)) cycle
      do iv = 1, nvp
        vmax = max(vmax, abs(dat(ip + (iv-1)*nph)))
      end do
    end do
    if (vmax <= 0) call fail('No nonzero unmasked line flux remains.')
    su = 0
    do ip = 1, nph
      if (.not. phase_valid(ip)) cycle
      do iv = 1, nvp
        su = su + (dat(ip + (iv-1)*nph)/vmax)**2
      end do
    end do
    flux_scale = vmax*sqrt(su/(count(phase_valid)*nvp))
    dat(1:nil) = dat(1:nil)/flux_scale
    eb(1:nil) = eb(1:nil)/flux_scale
    statistic_rank = 0
    do ip = 1, nph
      if (phase_valid(ip)) statistic_rank = statistic_rank + covariance_rank(ip)
      if (iw == 2) precision(:, :, ip) = precision(:, :, ip)*(flux_scale/normalization(ip))**2
    end do
    open(newunit=u, file='dop.normalization', status='replace')
    write(u, *) nph
    do ip = 1, nph
      write(u, *) normalization(ip)
    end do
    close(u)
  end subroutine statistical_input

  subroutine set_weights(eb, iw)
    ! Diagonal statistical weights, and the mean per-channel information
    ! (curvature of the likelihood), which sets the automatic alpha scale.
    real(dp), intent(in) :: eb(:)
    integer, intent(in) :: iw
    integer :: ip, iv, k, nused
    real(dp) :: total
    if (allocated(stat_weights)) deallocate(stat_weights)
    allocate(stat_weights(nil))
    stat_weights = 0
    total = 0
    nused = 0
    do ip = 1, nph
      if (.not. phase_valid(ip)) cycle
      do iv = 1, nvp
        k = ip + (iv-1)*nph
        stat_weights(k) = 1
        if (iw == 1) stat_weights(k) = 1/eb(k)**2
        ! Marginal precision 1/sigma^2, also for full covariance: the diagonal
        ! of an inverted near-singular covariance would overstate it.
        if (iw >= 1) then
          total = total + 1/eb(k)**2
        else
          total = total + 1
        end if
        nused = nused + 1
      end do
    end do
    information = likelihood_scale*total/nused
  end subroutine set_weights

  integer function edge_bin(x, edges)
    ! Largest i with edges(i) <= x (0 below the first edge).
    real(dp), intent(in) :: x, edges(:)
    integer :: lo, hi, mid
    lo = 0
    hi = size(edges)
    do while (hi - lo > 1)
      mid = (lo + hi)/2
      if (edges(mid) <= x) then
        lo = mid
      else
        hi = mid
      end if
    end do
    edge_bin = lo
  end function edge_bin

  real(dp) function normal_first(x, sigma)
    ! First antiderivative of the Gaussian CDF (of the step if sigma=0).
    real(dp), intent(in) :: x, sigma
    real(dp) :: z
    if (sigma <= 0) then
      normal_first = max(x, 0.0_dp)
    else
      z = x/sigma
      normal_first = x*0.5_dp*erfc(-z/sqrt(2.0_dp)) + sigma*exp(-0.5_dp*z*z)/sqrt(pi2)
    end if
  end function normal_first

  real(dp) function normal_second(x, sigma)
    ! Second antiderivative of the Gaussian CDF (of the step if sigma=0).
    real(dp), intent(in) :: x, sigma
    real(dp) :: z, cdf, pdf
    if (sigma <= 0) then
      normal_second = 0.5_dp*max(x, 0.0_dp)**2
    else
      z = x/sigma
      cdf = 0.5_dp*erfc(-z/sqrt(2.0_dp))
      pdf = exp(-0.5_dp*z*z)/sqrt(pi2)
      normal_second = 0.5_dp*((x*x + sigma*sigma)*cdf + x*sigma*pdf)
    end if
  end function normal_second

  real(dp) function square_cdf(x, ha, hb, sigma)
    ! CDF of a uniform square map pixel projected onto the line of sight:
    ! the convolution of boxes of half-widths ha and hb (a trapezoid),
    ! optionally convolved with a Gaussian line-spread function. Adjacent
    ! pixels tile exactly at every orbital phase.
    real(dp), intent(in) :: x, ha, hb, sigma
    real(dp) :: h1, h2, z
    h1 = max(ha, hb)
    h2 = min(ha, hb)
    if (h2 <= 1.0e-6_dp*h1) then
      z = (normal_first(x + h1, sigma) - normal_first(x - h1, sigma))/(2*h1)
    else
      z = (normal_second(x + h1 + h2, sigma) - normal_second(x + h1 - h2, sigma) &
           - normal_second(x - h1 + h2, sigma) + normal_second(x - h1 - h2, sigma)) &
          /(4*h1*h2)
    end if
    square_cdf = max(0.0_dp, min(1.0_dp, z))
  end function square_cdf

  subroutine ensure_capacity(needed)
    integer, intent(in) :: needed
    integer, allocatable :: index_tmp(:)
    real(dp), allocatable :: value_tmp(:)
    integer :: capacity
    if (needed <= size(proj_index)) return
    capacity = max(needed, 2*size(proj_index))
    allocate(index_tmp(capacity), value_tmp(capacity))
    index_tmp(1:size(proj_index)) = proj_index
    value_tmp(1:size(proj_value)) = proj_value
    call move_alloc(index_tmp, proj_index)
    call move_alloc(value_tmp, proj_value)
  end subroutine ensure_capacity

  subroutine build_projection(pha, vp)
    ! Exact square-pixel projection with exposure-phase quadrature.
    ! Phases are mid-exposure in radians; dpha holds full exposure widths.
    real(dp), intent(in) :: pha(:), vp(:)
    real(dp) :: vpl(nvp + 1), ww(nvp)
    real(dp) :: dvi, vx, vy, ph, sp, cp, vr, h1, h2, span, vfa, w
    integer :: ix, iy, ip, iq, iv, nq, ki, ka, klo, khi, ir, nnz, estimate
    va = max(abs(vp(1)), abs(vp(nvp)))
    if (nv < 4) call fail('Map size must be at least four pixels.')
    dvmap = 2*va/(nv - 1)
    dvi = minval(vp(2:nvp) - vp(1:nvp-1))
    if (dvi <= 0) call fail('Velocities must increase.')
    do iv = 1, nvp - 1
      vpl(iv + 1) = (vp(iv) + vp(iv + 1))/2
    end do
    vpl(1) = vp(1) - (vp(2) - vp(1))/2
    vpl(nvp + 1) = vp(nvp) + (vp(nvp) - vp(nvp - 1))/2
    if (allocated(indx)) deallocate(indx, indy, proj_ptr, proj_index, proj_value)
    allocate(indx(nv*nv), indy(nv*nv), proj_ptr(nv*nv + 1))
    estimate = max(1024, int(0.8_dp*nv*nv*nph*min(real(nvp, dp), 1.5_dp*dvmap/dvi + 3)))
    allocate(proj_index(estimate), proj_value(estimate))
    nnz = 0
    ir = 0
    do ix = 1, nv
      vx = dvmap*(ix - (nv + 1)/2.0_dp)
      do iy = 1, nv
        vy = dvmap*(iy - (nv + 1)/2.0_dp)
        if (vx**2 + vy**2 >= (va + dvmap/10)**2) cycle
        ir = ir + 1
        indx(ir) = ix
        indy(ir) = iy
        proj_ptr(ir) = nnz + 1
        do ip = 1, nph
          ww = 0
          klo = nvp + 1
          khi = 0
          ! Midpoint quadrature in phase, about three samples per bin of travel.
          nq = max(1, ceiling(abs(dpha(ip))*sqrt(vx*vx + vy*vy)/min(dvmap, dvi)*quadrature_density))
          nq = max(nq, phase_subsamples)
          if (nq > 100000) call fail('Exposure integration requires too many samples.')
          do iq = 1, nq
            ph = pha(ip) + dpha(ip)*((iq - 0.5_dp)/nq - 0.5_dp)
            sp = sin(ph)
            cp = cos(ph)
            vr = vx*sp - vy*cp
            h1 = 0.5_dp*dvmap*abs(sp)
            h2 = 0.5_dp*dvmap*abs(cp)
            span = h1 + h2 + 8*instrument_sigma
            ki = max(1, edge_bin(vr - span, vpl))
            ka = min(nvp, edge_bin(vr + span, vpl))
            vfa = 1 + af*exp(-min((vr/wid)**2, 50.0_dp))
            if (vfa <= 0) call fail('Absorption prescription makes the response nonpositive.')
            do iv = ki, ka
              w = square_cdf(vpl(iv + 1) - vr, h1, h2, instrument_sigma) &
                  - square_cdf(vpl(iv) - vr, h1, h2, instrument_sigma)
              ww(iv) = ww(iv) + max(w, 0.0_dp)/(nq*nph*vfa)
            end do
            klo = min(klo, ki)
            khi = max(khi, ka)
          end do
          if (khi < klo) cycle
          call ensure_capacity(nnz + khi - klo + 1)
          do iv = klo, khi
            if (ww(iv) <= 0) cycle
            nnz = nnz + 1
            proj_index(nnz) = ip + (iv - 1)*nph
            proj_value(nnz) = ww(iv)
          end do
        end do
      end do
    end do
    nol = ir
    proj_ptr(nol + 1) = nnz + 1
    if (allocated(dpx)) deallocate(dpx)
    allocate(dpx(nol, 4))
    dpx = 0
  end subroutine build_projection

  subroutine fproj(x, y)
    ! Forward projection: map (nol) -> data (nil).
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    integer :: j, k
    real(dp) :: xj
    y(1:nil) = 0
    do j = 1, nol
      xj = x(j)
      if (xj == 0) cycle
      do k = proj_ptr(j), proj_ptr(j + 1) - 1
        y(proj_index(k)) = y(proj_index(k)) + proj_value(k)*xj
      end do
    end do
  end subroutine fproj

  subroutine bproj(x, y)
    ! Adjoint (back) projection: data (nil) -> map (nol).
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    integer :: j, k
    real(dp) :: s
    do j = 1, nol
      s = 0
      do k = proj_ptr(j), proj_ptr(j + 1) - 1
        s = s + proj_value(k)*x(proj_index(k))
      end do
      y(j) = s
    end do
  end subroutine bproj

  subroutine baseline(dat, ih, datb, vmb)
    ! A flat circular map, scaled to baseline_fraction of the RMS flux at
    ! its projected peak, lets a positive latent map describe signed data.
    ! Gaussian objectives (ih=1) use this fixed level, so the returned map is
    ! bounded below by -vmb. The legacy log objective (ih=0) requires positive
    ! shifted data and raises the level when necessary (reported in status).
    real(dp), intent(in) :: dat(:)
    integer, intent(in) :: ih
    real(dp), intent(out) :: datb(:), vmb(:)
    real(dp) :: d0, bf
    integer :: i, ip
    vmb(1:nol) = 1
    call fproj(vmb, datb)
    d0 = maxval(datb(1:nil))
    if (d0 <= 0) call fail('The flat map has no projected flux.')
    bf = max(baseline_fraction, 1.0e-10_dp)
    baseline_raised = .false.
    if (ih == 0) then
      do i = 1, nil
        ip = mod(i - 1, nph) + 1
        if (.not. phase_valid(ip)) cycle
        if (datb(i) > 0) then
          if ((-dat(i) + 1.0e-10_dp)*d0/datb(i) > bf) then
            bf = (-dat(i) + 1.0e-10_dp)*d0/datb(i)
            baseline_raised = .true.
          end if
        end if
      end do
    end if
    datb(1:nil) = bf*datb(1:nil)/d0
    vmb(1:nol) = bf/d0
    baseline_level = bf
  end subroutine baseline

  subroutine gaussian_smear(x2, y2)
    ! Separable Gaussian exp(-r^2/ns^2), truncated at 1.5 ns, with
    ! half-sample reflection at the edges of the square array.
    real(dp), intent(in) :: x2(:,:)
    real(dp), intent(out) :: y2(:,:)
    real(dp) :: t(1 - max_half_width - 1:nv + max_half_width)
    integer :: i, j, h, n
    real(dp) :: s
    n = nv
    if (kernel_ns /= smoothing_width) then
      kernel_half = min(int(1.5*smoothing_width), max_half_width, n - 1)
      kernel_ns = smoothing_width
      kernel = 0
      s = 0
      do i = -kernel_half, kernel_half
        kernel(i) = exp(-real(i*i, dp)/real(smoothing_width, dp)**2)
        s = s + kernel(i)
      end do
      kernel(-kernel_half:kernel_half) = kernel(-kernel_half:kernel_half)/s
    end if
    h = kernel_half
    do i = 1, n
      do j = -h, 0
        t(j) = x2(i, 1 - j)
      end do
      t(1:n) = x2(i, 1:n)
      do j = n + 1, n + h
        t(j) = x2(i, 2*n - j + 1)
      end do
      do j = 1, n
        y2(i, j) = dot_product(t(j - h:j + h), kernel(-h:h))
      end do
    end do
    do j = 1, n
      do i = -h, 0
        t(i) = y2(1 - i, j)
      end do
      t(1:n) = y2(1:n, j)
      do i = n + 1, n + h
        t(i) = y2(2*n - i + 1, j)
      end do
      do i = 1, n
        y2(i, j) = dot_product(t(i - h:i + h), kernel(-h:h))
      end do
    end do
  end subroutine gaussian_smear

  subroutine chims(ps, chi)
    ! Default map: Gaussian-smoothed map, zero outside the circular support.
    real(dp), intent(in) :: ps(:)
    real(dp), intent(out) :: chi(:)
    real(dp) :: grid(nv, nv), smooth(nv, nv)
    integer :: j
    grid = 0
    do j = 1, nol
      grid(indx(j), indy(j)) = ps(j)
    end do
    call gaussian_smear(grid, smooth)
    do j = 1, nol
      chi(j) = smooth(indx(j), indy(j))
    end do
  end subroutine chims

  subroutine entr(ps, s, sp, chi, tem)
    ! S = -sum[x log(x/chi) - x + chi], chi = W x; sp = dS/dx including the
    ! dependence of chi on x (W is symmetric).
    real(dp), intent(in) :: ps(:)
    real(dp), intent(out) :: s, sp(:), chi(:), tem(:)
    real(dp) :: x(nol), y(nol)
    integer :: j
    call chims(ps, chi)
    s = 0
    do j = 1, nol
      tem(j) = log(ps(j)) - log(chi(j))
      s = s - ps(j)*(tem(j) - 1) - chi(j)
      sp(j) = -tem(j)
      x(j) = (ps(j) - chi(j))/chi(j)
    end do
    call chims(x, y)
    sp(1:nol) = sp(1:nol) + y
  end subroutine entr

  subroutine like(pht, ph, ih, h, dha)
    ! Likelihood H and its gradient with respect to the map.
    real(dp), intent(in) :: pht(:), ph(:)
    integer, intent(in) :: ih
    real(dp), intent(out) :: h, dha(:)
    real(dp) :: res(nil), phh(nil), pp, q
    integer :: i
    h = 0
    if (ih == 0 .and. stat_iw == 0) then
      do i = 1, nil
        phh(i) = 0
        if (stat_weights(i) == 0) cycle
        pp = max(ph(i), tiny(1.0_dp))
        h = h + pht(i)*log(pp) - pp
        phh(i) = pht(i)/pp - 1
      end do
    else
      res = pht(1:nil) - ph(1:nil)
      call apply_precision(res, phh, q)
      h = -0.5_dp*q
    end if
    h = h*likelihood_scale
    phh = phh*likelihood_scale
    call bproj(phh, dha)
  end subroutine like

  subroutine lamun(pht, ps, d1, d2, ph, chi, tem, al, ih, fla, fmu)
    ! Newton step in the plane of two directions, from exact first and
    ! second directional derivatives of Q = H + alpha S.
    real(dp), intent(in) :: pht(:), ps(:), d1(:), d2(:), ph(:), chi(:), tem(:), al
    integer, intent(in) :: ih
    real(dp), intent(out) :: fla, fmu
    real(dp) :: dph1(nil), dph2(nil), dch1(nol), dch2(nol)
    real(dp) :: rw(nil), z1(nil), z2(nil), res(nil)
    real(dp) :: hl, hll, hm, hlm, hmm, sl, sll, sm, slm, smm, q, q1, q2
    real(dp) :: pp, dpp, xl, xm, pl, pm, tl, tm, pz
    real(dp) :: ql, qll, qlm, qm, qmm, dd, denom, dr
    integer :: i, j
    call fproj(d1, dph1)
    call fproj(d2, dph2)
    hl = 0; hll = 0; hm = 0; hlm = 0; hmm = 0
    sl = 0; sll = 0; sm = 0; slm = 0; smm = 0
    if (ih == 0 .and. stat_iw == 0) then
      do i = 1, nil
        if (stat_weights(i) == 0) cycle
        pp = max(ph(i), tiny(1.0_dp))
        dpp = pht(i)/pp - 1
        xl = dph1(i)/pp
        xm = dph2(i)/pp
        hl = hl + dph1(i)*dpp
        hll = hll - xl**2*pht(i)
        hlm = hlm - xl*xm*pht(i)
        hm = hm + dph2(i)*dpp
        hmm = hmm - xm**2*pht(i)
      end do
    else
      res = pht(1:nil) - ph(1:nil)
      call apply_precision(res, rw, q)
      call apply_precision(dph1, z1, q1)
      call apply_precision(dph2, z2, q2)
      do i = 1, nil
        hl = hl + dph1(i)*rw(i)
        hm = hm + dph2(i)*rw(i)
        hlm = hlm - dph1(i)*z2(i)
      end do
      hll = -q1
      hmm = -q2
    end if
    hl = hl*likelihood_scale; hll = hll*likelihood_scale; hm = hm*likelihood_scale
    hlm = hlm*likelihood_scale; hmm = hmm*likelihood_scale
    if (al > 0) then
      call chims(d1, dch1)
      call chims(d2, dch2)
      do j = 1, nol
        pl = d1(j)/ps(j)
        pm = d2(j)/ps(j)
        tl = dch1(j)/chi(j)
        tm = dch2(j)/chi(j)
        pz = ps(j) - chi(j)
        sl = sl - d1(j)*tem(j) + pz*tl
        sm = sm - d2(j)*tem(j) + pz*tm
        sll = sll - ps(j)*(pl - tl)**2
        slm = slm - ps(j)*(pl - tl)*(pm - tm)
        smm = smm - ps(j)*(pm - tm)**2
      end do
    end if
    ql = hl + al*sl
    qll = hll + al*sll
    qlm = hlm + al*slm
    qm = hm + al*sm
    qmm = hmm + al*smm
    dd = qmm*qll - qlm**2
    denom = qll**2 + qmm**2
    fla = 0
    fmu = 1
    if (denom <= tiny(1.0_dp)) return
    dr = abs(dd)/denom
    if (dr < 1.0e-8_dp .or. dd <= 0 .or. qll >= 0) then
      if (qmm < -tiny(1.0_dp)) fmu = -qm/qmm
    else
      fmu = (qlm*ql - qll*qm)/dd
      fla = (qlm*qm - qmm*ql)/dd
    end if
  end subroutine lamun

  subroutine optimize(pht, al, ac, nim, ih, ps, nit, hs, s, ierr)
    ! Positivity-preserving preconditioned ascent of Q = H + alpha S, with
    ! two-direction Newton acceleration and an Armijo line search. Converged
    ! when ||x*grad|| / (||x*gH|| + ||x*alpha*gS||) <= ac.
    real(dp), intent(in) :: pht(:), al, ac
    integer, intent(in) :: nim, ih
    real(dp), intent(inout) :: ps(:)
    integer, intent(out) :: nit, ierr
    real(dp), intent(out) :: hs, s
    real(dp) :: ph(nil), dha(nol), dsa(nol), chi(nol), tem(nol)
    real(dp) :: ph_t(nil), dha_t(nol), dsa_t(nol), chi_t(nol), tem_t(nol)
    real(dp) :: dps(nol), dps0(nol), up(nol), trial(nol), grad(nol)
    real(dp) :: h, hnew, snew, qnew, su, sh, ss, denom, del, slope, step, fla, fmu
    integer :: j, ils
    logical :: accepted
    ierr = 1
    dps0 = 0
    dpx = 0
    s = 0
    dsa = 0
    call fproj(ps, ph)
    call like(pht, ph, ih, h, dha)
    if (al > 0) call entr(ps, s, dsa, chi, tem)
    do nit = 0, nim
      hs = h + al*s
      solver_objective = hs
      if (.not. ieee_is_finite(hs)) then
        ierr = 3
        return
      end if
      su = 0; sh = 0; ss = 0
      do j = 1, nol
        grad(j) = dha(j) + al*dsa(j)
        dps(j) = ps(j)*grad(j)
        su = su + dps(j)**2
        sh = sh + (ps(j)*dha(j))**2
        ss = ss + (ps(j)*al*dsa(j))**2
      end do
      denom = max(sqrt(sh) + sqrt(ss), 1.0e-12_dp)
      del = sqrt(su)/denom
      solver_stationarity = del
      if (nit == 0 .or. mod(nit, 10) == 0 .or. del <= ac) &
        write(log_unit, *) 'iteration objective stationarity ', nit, hs, del
      if (del <= ac .or. sqrt(su) <= 1.0e-12_dp) then
        ierr = 0
        return
      end if
      if (nit == nim) return
      if (nit == 0) dps0 = dps
      call lamun(pht, ps, dps0, dps, ph, chi, tem, al, ih, fla, fmu)
      up = fla*dps0 + fmu*dps
      slope = dot_product(grad, up)
      if (.not. ieee_is_finite(slope) .or. slope <= 0) then
        up = dps
        slope = dot_product(grad, up)
      end if
      step = 1
      do j = 1, nol
        if (up(j) < 0) step = min(step, -0.99_dp*ps(j)/up(j))
      end do
      accepted = .false.
      do ils = 1, 60
        trial = ps + step*up
        call fproj(trial, ph_t)
        call like(pht, ph_t, ih, hnew, dha_t)
        snew = 0
        dsa_t = 0
        if (al > 0) call entr(trial, snew, dsa_t, chi_t, tem_t)
        qnew = hnew + al*snew
        if (ieee_is_finite(qnew) .and. qnew >= hs + 1.0e-4_dp*step*slope) then
          accepted = .true.
          exit
        end if
        step = step/2
      end do
      if (.not. accepted) then
        ierr = 2
        return
      end if
      ! Reuse the accepted trial's projection, likelihood and entropy.
      ps = trial
      ph = ph_t
      dha = dha_t
      dsa = dsa_t
      h = hnew
      s = snew
      if (al > 0) then
        chi = chi_t
        tem = tem_t
      end if
      dps0 = step*up
      dpx(:, mod(nit, 4) + 1) = real(grad)
    end do
  end subroutine optimize

  real(dp) function rrms(dat, model)
    ! RMS statistic: sqrt(r^T W r / rank) over unmasked spectra.
    real(dp), intent(in) :: dat(:), model(:)
    real(dp) :: res(nil), weighted(nil), q
    res = dat(1:nil) - model(1:nil)
    call apply_precision(res, weighted, q)
    rrms = sqrt(max(q, 0.0_dp)/statistic_rank)
  end function rrms

  subroutine alpha_schedule(dat, al0, alf, nalpha, ac, nim, clim, ih, vmo, dato, &
                            alx, rr, sx, ier, nit_last, nit_total, steps, bracketed, plateau)
    ! Decrease alpha from al0 by factors alf until the RMS statistic meets
    ! clim. A bracketed crossing is refined by bisection in log(alpha), so
    ! the accepted map is close to the most regularized one meeting clim.
    ! Stops early (ier=1, plateau=1) when the statistic no longer improves.
    ! ier: 0 success; 1 target not met; 10+k optimizer failure k.
    real(dp), intent(in) :: dat(:), al0, alf, ac, clim
    integer, intent(in) :: nalpha, nim, ih
    real(dp), intent(inout) :: vmo(:)
    real(dp), intent(out) :: dato(:), alx, rr, sx
    integer, intent(out) :: ier, nit_last, nit_total, steps, bracketed, plateau
    real(dp) :: vm_try(nol), vm_hi(nol), model(nil)
    real(dp) :: al, al_hi, al_lo, al_mid, rr_try, rr_prev, drop, drop_prev, ratio, ratio_prev, hs, s
    integer :: ial, ib, nit, ierr
    ier = 1
    bracketed = 0
    plateau = 0
    nit_total = 0
    nit_last = 0
    steps = 0
    rr = huge(1.0_dp)
    rr_prev = huge(1.0_dp)
    drop_prev = -1
    ratio_prev = -1
    rms_floor_estimate = -1
    alx = al0
    sx = 0
    vm_try = vmo
    al_hi = al0
    do ial = 1, nalpha
      al = al0/alf**(ial - 1)
      call solve(al, ierr)
      if (ierr /= 0) then
        ! Report the failing solve; no map is accepted from it.
        rr = rr_try
        alx = al
        nit_last = nit
        return
      end if
      if (rr_try <= clim) then
        call accept(al)
        ier = 0
        if (ial > 1) then
          bracketed = 1
          al_lo = al
          do ib = 1, 8
            if (clim - rr <= 0.005_dp*clim) exit
            al_mid = sqrt(al_hi*al_lo)
            vm_try = vm_hi
            call solve(al_mid, ierr)
            if (ierr /= 0) then
              ier = 0
              exit
            end if
            if (rr_try <= clim) then
              call accept(al_mid)
              al_lo = al_mid
            else
              vm_hi = vm_try
              al_hi = al_mid
            end if
          end do
        end if
        return
      end if
      if (ial > 1) then
        drop = rr_prev - rr_try
        ! Stop when the statistic has stopped improving, or when geometric
        ! extrapolation of steadily shrinking improvements (a lower bound in
        ! practice) shows the remaining schedule cannot reach the target.
        if (drop < 1.0e-3_dp*rr_prev) then
          rms_floor_estimate = rr_try
          call accept(al)
          plateau = 1
          ier = 1
          return
        end if
        if (drop_prev > 0) then
          ratio = drop/drop_prev
          ! Two consecutive shrinking steps are required, and the slower of
          ! the two decay rates is used, so the extrapolated floor errs low.
          if (ratio > 0 .and. ratio < 0.8_dp .and. ratio_prev > 0 .and. ratio_prev < 0.8_dp) then
            rms_floor_estimate = rr_try - drop*max(ratio, ratio_prev)/(1 - max(ratio, ratio_prev))
            if (rms_floor_estimate > clim*(1 + 1.0e-3_dp)) then
              call accept(al)
              plateau = 1
              ier = 1
              return
            end if
          end if
          ratio_prev = ratio
        end if
        drop_prev = drop
      end if
      call accept(al)
      ier = 1
      vm_hi = vm_try
      al_hi = al
      rr_prev = rr_try
    end do

  contains

    subroutine solve(alpha, status)
      real(dp), intent(in) :: alpha
      integer, intent(out) :: status
      call optimize(dat, alpha, ac, nim, ih, vm_try, nit, hs, s, status)
      steps = steps + 1
      nit_total = nit_total + nit
      call fproj(vm_try, model)
      rr_try = rrms(dat, model)
      write(log_unit, *) 'ni, al, hs, rr:', nit, alpha, hs, rr_try, clim
      if (status /= 0 .or. .not. ieee_is_finite(rr_try)) then
        ier = 10 + status
        if (status == 0) ier = 19
        write(log_unit, *) 'optimizer failure ', ier
        status = max(status, 1)
      end if
    end subroutine solve

    subroutine accept(alpha)
      real(dp), intent(in) :: alpha
      vmo(1:nol) = vm_try
      dato(1:nil) = model
      alx = alpha
      rr = rr_try
      sx = s
      nit_last = nit
    end subroutine accept

  end subroutine alpha_schedule

end module dop_science
