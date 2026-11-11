! Doppler tomography driver: reads dop.in and dopin, writes dop.out,
! dop.status, dop.log and dop.normalization. See dop_core.f90.
program dopmap
  use dop_science
  implicit none
  integer :: ih, iw, ns, nim, nal, ipri, norm, map_size, iph, i, ip, iv, j, m, u
  integer :: ier, nit, nit_total, steps, bracketed, plateau, iconverged, nvalid
  real(dp) :: pb0, pb1, ac, al0, alf, clim, w0, al, rr, sx, sdat, s1, s2, fident
  real(dp) :: dev, denu, alpha_start, t1, t2
  character(len=256) :: idat
  real(dp), allocatable :: pha(:), vp(:), dat(:), eb(:), datb(:), vmb(:), vmo(:), dato(:)
  real(dp), allocatable :: reco(:,:), dpxu(:,:)

  open(newunit=log_unit, file='dop.log', status='replace')
  call read_parameters()
  smoothing_width = ns
  nv = map_size

  open(newunit=u, file='dopin', status='old', action='read')
  read(u, *) nph, nvp, w0
  read(u, '(a)') idat
  nil = nph*nvp
  allocate(pha(nph), dpha(nph), vp(nvp), dat(nil), eb(nil))
  read(u, *) (pha(i), i = 1, nph)
  read(u, *) iph
  if (iph == 1) then
    read(u, *) (dpha(i), i = 1, nph)
  else
    dpha = 1.0_dp/nph
  end if
  pha = pha*pi2
  dpha = dpha*pi2
  read(u, *) (vp(i), i = 1, nvp)
  read(u, *) (dat(i), i = 1, nil)
  eb = 1
  if (iw >= 1) read(u, *) (eb(i), i = 1, nil)
  close(u)

  stat_iw = iw
  call set_phase_mask(pha, pb0*pi2, pb1*pi2)
  call statistical_input(dat, eb, iw, norm)
  if (iw >= 1) ih = 1
  ! Legacy (unweighted) objectives are scaled per spectral channel, so alpha
  ! does not depend on spectral sampling and keeps its historical meaning.
  likelihood_scale = 1
  if (iw == 0) likelihood_scale = 1.0_dp/nvp

  call cpu_time(t1)
  call build_projection(pha, vp)
  call cpu_time(t2)
  write(log_unit, '(a,f10.2,a,i12)') ' cpu for geometry', t2 - t1, '   nonzeros', proj_ptr(nol + 1) - 1

  allocate(datb(nil), vmb(nol), vmo(nol), dato(nil))
  call baseline(dat, ih, datb, vmb)
  dat = dat + datb
  sdat = 0
  nvalid = 0
  do ip = 1, nph
    if (.not. phase_valid(ip)) cycle
    nvalid = nvalid + 1
    do iv = 1, nvp
      sdat = sdat + dat(ip + (iv - 1)*nph)
    end do
  end do
  sdat = sdat*nph/nvalid
  call set_weights(eb, iw)

  ! Automatic start: the legacy alpha 0.002 per unit channel information.
  alpha_start = al0
  if (al0 <= 0) alpha_start = 0.002_dp*nvp*information
  write(log_unit, *) 'alpha start, information, baseline ', alpha_start, information, baseline_level

  ! Flat initial map with the data's total flux.
  do j = 1, nol
    vmo(j) = max(sdat/nol, vmb(j))
  end do
  call alpha_schedule(dat, alpha_start, alf, max(1, nal), ac, nim, clim, ih, vmo, dato, &
                      al, rr, sx, ier, nit, nit_total, steps, bracketed, plateau)

  iconverged = 0
  if (ier == 0 .and. ieee_is_finite(rr)) iconverged = 1
  open(newunit=u, file='dop.status', status='replace')
  write(u, *) 'converged ', iconverged
  write(u, *) 'optimizer_status ', ier
  write(u, *) 'fit_statistic ', rr
  write(u, *) 'fit_target ', clim
  write(u, *) 'iterations ', nit
  write(u, *) 'total_iterations ', nit_total
  write(u, *) 'alpha ', al
  write(u, *) 'alpha_start ', alpha_start
  write(u, *) 'alpha_steps ', steps
  write(u, *) 'target_bracketed ', bracketed
  write(u, *) 'rms_plateau ', plateau
  write(u, *) 'rms_floor_estimate ', rms_floor_estimate
  write(u, *) 'entropy ', sx
  write(u, *) 'objective ', solver_objective
  write(u, *) 'stationarity ', solver_stationarity
  write(u, *) 'flux_scale ', flux_scale
  write(u, *) 'likelihood ', ih
  write(u, *) 'weight_mode ', iw
  write(u, *) 'statistic_rank ', statistic_rank
  write(u, *) 'information ', information
  write(u, *) 'baseline ', baseline_level
  write(u, *) 'baseline_raised ', merge(1, 0, baseline_raised)
  close(u)
  if (iconverged /= 1) stop 1
  if (.not. all(ieee_is_finite(vmo(1:nol)))) stop 1

  ! Identification number (historical).
  s1 = 0
  s2 = 0
  m = max(1, nol/30)
  do i = 1, nol, m
    s1 = s1 + vmo(i)
  end do
  m = max(1, nil/30)
  do i = 1, nil, m
    s2 = s2 + dat(i)
  end do
  fident = 0
  if (abs(s2) > tiny(1.0_dp)) fident = sin(100*(s1/s2))

  ! Map density in input spectral-density units per km/s: its velocity-area
  ! integral equals the line integral over radial velocity in km/s.
  dev = (vp(nvp) - vp(1))/(nvp - 1)
  denu = dev*1.0e5_dp
  do j = 1, nol
    vmo(j) = (vmo(j) - vmb(j))*flux_scale/dvmap**2*denu/nph
  end do
  ! Restore spectra to their original physical flux scale for diagnostics.
  do ip = 1, nph
    do iv = 1, nvp
      i = ip + (iv - 1)*nph
      dat(i) = (dat(i) - datb(i))*flux_scale/normalization(ip)
      dato(i) = (dato(i) - datb(i))*flux_scale/normalization(ip)
    end do
  end do
  allocate(reco(nv, nv), dpxu(nv, nv))
  reco = 0
  dpxu = 0
  do j = 1, nol
    reco(indx(j), indy(j)) = vmo(j)
  end do

  open(newunit=u, file='dop.out', status='replace')
  write(u, *) nph, nvp, nv, w0, fident
  write(u, *) trim(idat)
  write(u, '(1p8e24.15E3)') (pha(i), i = 1, nph)
  write(u, *) iph
  if (iph == 1) write(u, '(1p8e24.15E3)') (dpha(i), i = 1, nph)
  write(u, '(1p8e24.15E3)') (vp(i), i = 1, nvp)
  write(u, '(1p8e24.15E3)') (dat(i), i = 1, nil)
  write(u, *) ih, iw, pb0, pb1, ns, ac, al, clim, norm, wid, af
  write(u, *) nv, va, ' dopmap '
  write(u, '(1p6e24.15E3)') ((reco(i, j), j = 1, nv), i = 1, nv)
  write(u, '(i12,a)') nil, ' reconstructed data '
  write(u, '(1p6e24.15E3)') (dato(i), i = 1, nil)
  write(u, '(i12,a)') nv, ' 4  change maps '
  do m = 1, 4
    do j = 1, nol
      dpxu(indx(j), indy(j)) = dpx(j, m)
    end do
    write(u, '(1p6e24.15E3)') ((dpxu(i, j), j = 1, nv), i = 1, nv)
  end do
  close(u)
  close(log_unit)

contains

  subroutine read_parameters()
    integer :: pu, ios
    open(newunit=pu, file='dop.in', status='old', action='read', iostat=ios)
    if (ios /= 0) call fail('Cannot open dop.in.')
    read(pu, *) ih
    read(pu, *) iw
    read(pu, *) pb0, pb1
    read(pu, *) ns
    read(pu, *) ac
    read(pu, *) nim
    read(pu, *) al0, alf, nal
    read(pu, *) clim
    read(pu, *) ipri
    read(pu, *) norm
    read(pu, *) wid, af
    read(pu, *) instrument_sigma, phase_subsamples, baseline_fraction
    read(pu, *, iostat=ios) map_size
    if (ios /= 0) call fail('dop.in lacks the map size; refresh the Fortran sources.')
    close(pu)
    if (pb1 < pb0 .or. pb1 - pb0 >= 1) call fail('Require 0 <= pb1 - pb0 < 1.')
    if (ih < 0 .or. ih > 1 .or. iw < 0 .or. iw > 2) call fail('Invalid ih or iw.')
    if (ns < 1 .or. ac <= 0 .or. nim < 1 .or. alf <= 1 .or. clim <= 0) &
      call fail('Invalid ns, ac, nim, alf or clim.')
    if (norm < 0 .or. norm > 1) call fail('norm must be 0 or 1.')
    if (instrument_sigma < 0 .or. phase_subsamples < 0 .or. baseline_fraction < 0) &
      call fail('Invalid LSF sigma, phase subsamples or baseline fraction.')
    if (map_size < 4 .or. mod(map_size, 2) /= 0) call fail('Map size must be even and >= 4.')
    if (wid <= 0) then
      wid = 1
      af = 0
    end if
    write(log_unit, '(a)') ' Floating-default maximum entropy Doppler tomography'
    write(log_unit, '(a,2i3)') ' ih, iw    ', ih, iw
    write(log_unit, '(a,2f8.4)') ' pb0, pb1  ', pb0, pb1
    write(log_unit, '(a,i5)') ' ns        ', ns
    write(log_unit, '(a,1pe10.2)') ' ac        ', ac
    write(log_unit, '(a,i8)') ' nim       ', nim
    write(log_unit, '(a,1p2e12.4,i5)') ' al0 alf nal', al0, alf, nal
    write(log_unit, '(a,1pe12.4)') ' clim      ', clim
    write(log_unit, '(a,i3)') ' norm      ', norm
    write(log_unit, '(a,1p2e12.4)') ' wid, af   ', wid, af
    write(log_unit, '(a,i6)') ' map size  ', map_size
  end subroutine read_parameters

end program dopmap
