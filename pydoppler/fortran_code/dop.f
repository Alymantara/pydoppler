      module dop_science
      implicit none
      integer :: stat_iw=0, nphase=0, nvelocity=0
      integer :: phase_subsamples=0, statistic_rank=0
      real*8 :: flux_scale=1.d0, instrument_sigma=0.d0
      real*8 :: baseline_fraction=0.1d0
      real*8 :: solver_objective=0.d0,solver_stationarity=0.d0
      real*8, allocatable :: precision(:,:,:),stat_weights(:)
      real*8, allocatable :: normalization(:)
      integer, allocatable :: covariance_rank(:)
      logical, allocatable :: phase_valid(:)
      contains
      subroutine apply_precision(residual,weighted,quadratic)
      implicit none
      real*8, intent(in) :: residual(*)
      real*8, intent(out) :: weighted(*),quadratic
      real*8 :: rv(nvelocity),zv(nvelocity)
      integer :: ip,iv,k
      quadratic=0.d0
      do ip=1,nphase
        do iv=1,nvelocity
          k=ip+(iv-1)*nphase
          rv(iv)=residual(k)
          weighted(k)=0.d0
        enddo
        if (.not.phase_valid(ip)) cycle
        if (stat_iw.eq.2) then
          zv=matmul(precision(:,:,ip),rv)
        else
          do iv=1,nvelocity
            k=ip+(iv-1)*nphase
            zv(iv)=stat_weights(k)*rv(iv)
          enddo
        endif
        do iv=1,nvelocity
          k=ip+(iv-1)*nphase
          weighted(k)=zv(iv)
        enddo
        quadratic=quadratic+dot_product(rv,zv)
      enddo
      end subroutine
      end module dop_science

      program dopmap
      use dop_science
      use, intrinsic :: ieee_arithmetic
c Based on H. C. Spruit's 1999 floating-default maximum-entropy code.
c Revised likelihood, optimization, quadrature and physical output units.
      implicit real*8 (a-h,o-z)
      real*4 dpx
      include 'emap.par'
      common/phases/dpha(npm)
c the 'd' lines are the nonsense IBM machines need to make them stop on an overflow.
c      include '/usr/include/fpdc.h'
c      include '/usr/include/fexcp.h'
      dimension dat(nd),pha(npm),vp(nvpm)
      real*4 p(ndi,nr),pt(nri,nd)
      dimension vmo(nr),vmb(nr),dato(nd),datb(nd),eb(nd)
      dimension reco(nvm,nvm),dpxu(nvm,nvm)
      character*20 spec
      character*80 idat
      common /folds/ nph,nvp
      common /foldn/ nv,indx(nr),indy(nr),nol
      common /difpix/ dpx(nr,4)
      common /cabsf/ wid,af
c      fpstat(fpve)=.true.
c      fpstat(fpoe)=.true.
c      fpstat(fpue)=.true.
c      fpstat(fpze)=.true.
c      fpstat(fpxe)=.false.
c      call fpsets(fpstat)
c      call signal(sigtrap,xl__trce)
      pi2=8.d0*atan(1.d0)

c read parameter file
      call repar(ns,ac,nim,
     *     al0,alf,nal,clim,ipri,spec,ih,iw,pb0,pb1,norm)
c read input data
      open(4,file=spec)
      read(4,*) nph,nvp,w0
      read(4,'(a)') idat
      if (nph.ne.npm) call erri('nph.ne.npm',nph,npm)
      if (nvp.ne.nvpm) call erri('nvp.ne.npvm',nvp,nvpm)
      nv=nvm
      nil=nvp*nph
      read(4,*) (pha(i),i=1,nph)
      read(4,*) iph
      if (iph.eq.1) then
         read(4,*) (dpha(i),i=1,nph)
      else
         do i=1,nph
           dpha(i)=1./nph
         enddo
      endif
      do i=1,nph
        pha(i)=pha(i)*pi2
        dpha(i)=dpha(i)*pi2
      enddo
c make blanking array corresponding to (pb0,pb1)
      call blanks(pha,nph,pb0,pb1)
      read(4,*) (vp(i),i=1,nvp)
      read(4,*) (dat(i),i=1,nil)
      if (iw.ge.1) read(4,*) (eb(i),i=1,nil)
      close(4)
c Set one explicitly shared statistical model for likelihood and diagnostics.
      stat_iw=iw
      nphase=nph
      nvelocity=nvp
      allocate(stat_weights(nil),normalization(nph))
      allocate(phase_valid(nph),covariance_rank(nph))
      call statistical_input(dat,eb,iw,norm,nph,nvp)
      if (iw.ge.1) ih=1
      call cpu_time(t1)
      call prom(va,nv,pha,nph,vp,nvp,pt,p)
      call cpu_time(t2)
      write(10,'('' cpu for geometry'',f10.2)')t2-t1
      if (ipri.gt.0) write(*,'('' cpu for geometry'',f10.2)')t2-t1
      t1=t2
c A deterministic projected baseline permits signed continuum residuals.
      call basel(dat,nph,nvp,baseline_fraction,bf,nv,va,
     *           pt,datb,vmb)
      do i=1,nil
        dat(i)=dat(i)+datb(i)
      enddo
      sdat=0
      nvalid=0
      do ip=1,nph
        if (.not.phase_valid(ip)) cycle
        nvalid=nvalid+1
        do iv=1,nvp
          sdat=sdat+dat(ip+(iv-1)*nph)
        enddo
      enddo
      sdat=sdat*nph/nvalid
      call weighd(dat,datb,pha,nph,nvp,eb,iw)

c reconstruction
c initialize (flat)
      do j=1,nol
        vmo(j)=sdat/nol
      enddo
c decreasing sequence of alfas until required fit to input data
      call alfs(dat,eb,iw,p,pt,al0,alf,nal,ac,nim,vmo,
     *           dato,nil,nol,nit,al,rr,ns,ih,clim,ipri,ier,t2,sx)
      open(12,file='dop.status',status='replace')
      iconverged=0
      if (ier.eq.0.and.ieee_is_finite(rr)) iconverged=1
      write(12,*) 'converged ',iconverged
      write(12,*) 'optimizer_status ',ier
      write(12,*) 'fit_statistic ',rr
      write(12,*) 'fit_target ',clim
      write(12,*) 'iterations ',nit
      write(12,*) 'alpha ',al
      write(12,*) 'entropy ',sx
      write(12,*) 'objective ',solver_objective
      write(12,*) 'stationarity ',solver_stationarity
      write(12,*) 'flux_scale ',flux_scale
      write(12,*) 'likelihood ',ih
      write(12,*) 'weight_mode ',iw
      write(12,*) 'statistic_rank ',statistic_rank
      close(12)
      if (iconverged.ne.1) stop 1
      if (.not.all(ieee_is_finite(vmo(1:nol)))) stop 1
c make an identification number
      s1=0
      s2=0
      m=max(1,nol/30)
      do i=1,nol,m
        s1=s1+vmo(i)
      enddo
      m=max(1,nil/30)
      do i=1,nil,m
        s2=s2+dat(i)
      enddo
      fident=0.d0
      if (abs(s2).gt.tiny(1.d0)) fident=sin(100*(s1/s2))
c Map density: input spectral-density units per km/s.
c Its velocity-area integral equals the line integral over radial km/s.
c velocity grid size
      dv=2*va/(nv-1)
      delv2=dv**2
c Spectral velocity interval (stored in cm/s).
      dev=(vp(nvp)-vp(1))/(nvp-1)
      denu=dev*1.d5
      do i=1,nol
        vmo(i)=(vmo(i)-vmb(i))*flux_scale/delv2*denu/nph
      enddo
c Restore spectra to their original physical flux scale for diagnostics.
      do ip=1,nph
        do iv=1,nvp
          i=ip+(iv-1)*nph
          dat(i)=(dat(i)-datb(i))*flux_scale/normalization(ip)
          dato(i)=(dato(i)-datb(i))*flux_scale/normalization(ip)
          datb(i)=0.d0
        enddo
      enddo
c make unfolded reconstructed image
      reco=0.d0
      dpxu=0.d0
      do j=1,nol
        reco(indx(j),indy(j))=vmo(j)
      enddo

      open(3,file='dop.out')
c write array sizes, w0, identification
      write(3,*) nph,nvp,nvm,w0,fident
c copy second line of input data file
      write(3,*)idat
c copy phase and velocity coordinates of input spectrum
      write(3,'(1p8e24.15E3)')(pha(i),i=1,nph)
      write(3,*)iph
      if (iph.eq.1) write(3,'(1p8e24.15E3)')(dpha(i),i=1,nph)
      write(3,'(1p8e24.15E3)')(vp(i),i=1,nvp)
c write input spectrum
      write(3,'(1p8e24.15E3)')(dat(i)-datb(i),i=1,nil)
      write(3,*)
     *      ih,iw,pb0/pi2,pb1/pi2,ns,ac,al,clim,norm,wid,af
c write dopmap and reco
      write(3,*) nv,va,' dopmap '
      write(3,'(1p6e24.15E3)')((reco(i,j),j=1,nv),i=1,nv)
      write(3,'(i6,'' reconstructed data '')')nil
      write(3,'(1p6e24.15E3)')(dato(i)-datb(i),i=1,nil)
      write(3,'(i6,'' 4  change maps '')')nv
      do ii=1,4
        do j=1,nol
          dpxu(indx(j),indy(j))=dpx(j,ii)
        enddo
        write(3,'(1p6e24.15E3)')((dpxu(i,j),j=1,nv),i=1,nv)
      enddo
      close(3)
      end

c      block data
c      include '/usr/include/fpdc.h'
c      include '/usr/include/fpdt.h'
c      end

      subroutine statistical_input(dat,eb,iw,norm,nph,nvp)
      use dop_science
      use, intrinsic :: ieee_arithmetic
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension dat(*),eb(*)
      common/blank/bl(npm)
      logical bl
      phase_valid=.not.bl
      normalization=1.d0
      covariance_rank=nvp
      if (count(phase_valid).eq.0) then
        write(*,*) 'No unmasked orbital phases remain.'
        stop 1
      endif
      if (iw.lt.0.or.iw.gt.2) stop 1
      do i=1,nph*nvp
        if (.not.ieee_is_finite(dat(i))) stop 1
        if (iw.ge.1) then
          if (.not.ieee_is_finite(eb(i)).or.eb(i).le.0) stop 1
        else
          eb(i)=1.d0
        endif
      enddo
      if (iw.eq.2) then
        allocate(precision(nvp,nvp,nph))
        open(13,file='dop.precision',status='old',iostat=ios)
        if (ios.ne.0) stop 1
        read(13,*,iostat=ios) npcheck,nvcheck
        if (ios.ne.0.or.npcheck.ne.nph.or.nvcheck.ne.nvp) stop 1
        do ip=1,nph
          read(13,*,iostat=ios) covariance_rank(ip)
          if (ios.ne.0.or.covariance_rank(ip).lt.1.or.
     *        covariance_rank(ip).gt.nvp) stop 1
          read(13,*,iostat=ios)
     *       ((precision(i,j,ip),j=1,nvp),i=1,nvp)
          if (ios.ne.0) stop 1
          if (.not.all(ieee_is_finite(precision(:,:,ip)))) stop 1
        enddo
        close(13)
      endif
c Optional exact normalization uses the continuum-subtracted line integral.
      if (norm.eq.1) then
        avg=0.d0
        do ip=1,nph
          if (.not.phase_valid(ip)) cycle
          su=0.d0
          do iv=1,nvp
            su=su+dat(ip+(iv-1)*nph)
          enddo
          if (su.le.0.d0) then
            write(*,*) 'Phase normalization requires positive line flux.'
            stop 1
          endif
          normalization(ip)=su
          avg=avg+su
        enddo
        avg=avg/count(phase_valid)
        do ip=1,nph
          if (phase_valid(ip)) normalization(ip)=avg/normalization(ip)
        enddo
      endif
      do ip=1,nph
        do iv=1,nvp
          k=ip+(iv-1)*nph
          dat(k)=dat(k)*normalization(ip)
          eb(k)=eb(k)*normalization(ip)
        enddo
      enddo
c Stable RMS scale makes every subsequent numerical choice unit invariant.
      vmax=0.d0
      do ip=1,nph
        if (.not.phase_valid(ip)) cycle
        do iv=1,nvp
          vmax=max(vmax,abs(dat(ip+(iv-1)*nph)))
        enddo
      enddo
      if (vmax.le.0.d0) then
        write(*,*) 'No nonzero unmasked line flux remains.'
        stop 1
      endif
      su=0.d0
      do ip=1,nph
        if (.not.phase_valid(ip)) cycle
        do iv=1,nvp
          su=su+(dat(ip+(iv-1)*nph)/vmax)**2
        enddo
      enddo
      flux_scale=vmax*sqrt(su/(count(phase_valid)*nvp))
      do i=1,nph*nvp
        dat(i)=dat(i)/flux_scale
        eb(i)=eb(i)/flux_scale
      enddo
      statistic_rank=0
      do ip=1,nph
        if (phase_valid(ip)) statistic_rank=statistic_rank+
     *      covariance_rank(ip)
        if (iw.eq.2) precision(:,:,ip)=precision(:,:,ip)*
     *      (flux_scale/normalization(ip))**2
      enddo
      open(13,file='dop.normalization',status='replace')
      write(13,*) nph
      do ip=1,nph
        write(13,*) normalization(ip)
      enddo
      close(13)
      end

      subroutine basel(dat,nph,nvp,bfa,bf,nv,va,pt,datb,vmb)
c compute a baseline flux to get a flat background level in dopmap
c datb is the data corresponding to a flat (circular) dopmap
c vmb is the corresponding flat dopmap
      use dop_science
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      common /foldn/ nvf,indx(nr),indy(nr),nol
      dimension dat(*),datb(*),vmb(*)
      real*4 pt(nri,nd)
      if (nol.gt.nr) call erri('nol>nr',nol,nr)

c make circle of radius va on a zero background
      fv2=nv/2.
      dv=2*va/(nv-1)
      do ir=1,nol
        ix=indx(ir)
        iy=indy(ir)
        vx=dv*(ix-fv2-0.5)
        vy=dv*(iy-fv2-0.5)
        v=sqrt(vx**2+vy**2)
        vmb(ir)=1.d0
      enddo

      nil=nph*nvp
c produce spectrum from this map
      call fproj(pt,nil,nol,vmb,datb)
c Normalize the projected flat image by its largest positive bin.
      d0=maxval(datb(1:nil))
      if (d0.le.0.d0) stop 1
      bf=bfa
      do i=1,nil
        ip=mod(i-1,nph)+1
        if (.not.phase_valid(ip)) cycle
        if (datb(i).gt.0.d0) then
          bf=max(bf,(-dat(i)+1.d-10)*d0/datb(i))
        else if (dat(i).lt.0.d0) then
          stop 1
        endif
      enddo
      bf=max(bf,1.d-10)
      do i=1,nil
        datb(i)=bf*datb(i)/d0
      enddo
      do j=1,nol
        vmb(j)=vmb(j)*bf/d0
      enddo
      return
      end

      subroutine blanks(pha,nph,pb0,pb1)
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension pha(*)
      common/blank/bl(npm)
      logical bl
      pi2=8.d0*atan(1.d0)
      width=(pb1-pb0)/pi2
      tolerance=16.d0*epsilon(1.d0)
      do i=1,nph
        distance=modulo((pha(i)-pb0)/pi2,1.d0)
        bl(i)=distance.gt.tolerance.and.distance.lt.width-tolerance
      enddo
      end

      subroutine alfs(dat,eb,iw,p,pt,al0,alf,nal,ac,nim,vmo,
     *           dato,nil,nol,nit,alx,rr,ns,ih,clim,ipri,ier,t2,sx)
      use, intrinsic :: ieee_arithmetic
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension vm1(nr),vmo(*),dat(nd),eb(nd),dat1(nd),dato(nd)
      real*4 p(nfp,nr),pt(nrp,nf)
      vm1(1:nol)=vmo(1:nol)
      ier=1
      nit=0
      rr=huge(1.d0)
      sx=0.d0
      alx=al0
      nal1=max(1,abs(nal))
      call cpu_time(t1)
      do ial=1,nal1
        al=al0/alf**(ial-1)
        aci=ac
        call rldop(dat,p,pt,al,aci,1.d5,nim,vm1,nil,nol,
     *             nit,hs,ns,ih,ipri,sf,ierr)
        call fproj(pt,nil,nol,vm1,dat1)
        rr=rrms(dat,eb,iw,dat1,nil)
        alx=al
        sx=sf
        write(10,*) 'ni, al, hs, rr:',nit,al,hs,rr,clim
        call cpu_time(t2)
        write(10,*) 'cpu for iteration ',t2-t1
        t1=t2
        if (ierr.ne.0.or..not.ieee_is_finite(rr)) then
          ier=10+ierr
          if (ierr.eq.0) ier=19
          write(10,*) 'optimizer failure ',ier
          return
        endif
c Only an optimizer-converged solution can replace the previous model.
        vmo(1:nol)=vm1(1:nol)
        dato(1:nil)=dat1(1:nil)
        if (rr.le.clim) then
          ier=0
          return
        endif
      enddo
      end

      subroutine repar(ns,ac,nim,
     *     al0,alf,nal,clim,ipri,spec,ih,iw,pb0,pb1,norm)
c read input parameters
      use dop_science
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      character*(*) spec
      character*40 txt
c output channels
      dimension kch(2)
      common /fold/ nvx,nvy
      common /cabsf/ wid,af
      pi2=8.d0*atan(1.d0)
c input data file
      open(4,file='dop.in')
c type of likelihood to be used (0=log, 1=weighted log, 2=rms, 3=weighted rms)
      read(4,*)ih
c read weights from lightcurve file if iw=1
      read(4,*)iw
c read range of phases to be ignored in mapping and convert to radians
      read(4,*)pb0,pb1
      if(pb1-1.gt.pb0) call errb('all phases ignored:',pb0,pb1)
      if(pb1.lt.pb0) call errb('pb0 must be <pb1',pb0,pb1)
      na=nvm
      if (mod(na,2).ne.0) call erri('nvm not even',na,0)
      nvx=na
      nvy=na
c smearing width in default map
      read(4,*)ns
      if (ns.le.0) ns=4
c accuracy of RL iteration
      read(4,*)ac
      if (ac.le.0.) ac=3e-3
c max no of iterations
      read(4,*)nim
      if (nim.le.0) nim=100
c first value of alpha parameter for penalty function
      read(4,*)al0,alf,nal
      if (al0.le.0.)
     *     al0=0.003*(float(nvm)/nvpm)**2*min(npm,nvpm)/float(nvpm)
      if (alf.le.0.) alf=1.7
c number of alphas to be tried (if <0, only first value is used)
      if (nal.eq.0) nal=-1
      read(4,*)clim
      if (clim.le.0.) clim=1.6
c printout control
      read(4,*)ipri
      if (ipri.lt.1) ipri=1
      if (ipri.gt.2) ipri=2
c read name of input spectrum
      spec='dopin'
c spectra to be normalized to constant wavelength-integrated light?
      read(4,*)norm
c central absorption parameters
      read(4,*)wid,af
      if (wid.eq.0.) wid=2.5e7
      if (wid.le.0.) then
        wid=1
        af=0
      endif
      read(4,*,iostat=ios) instrument_sigma,phase_subsamples,
     *                     baseline_fraction
      if (ios.ne.0) then
        instrument_sigma=0.d0
        phase_subsamples=0
        baseline_fraction=0.1d0
      endif
      if (instrument_sigma.lt.0.or.phase_subsamples.lt.0.or.
     *    baseline_fraction.lt.0) stop 1
      close (4)
      kch(2)=6
      kch(1)=10
      open(kch(1),file='dop.log')
      if (ipri.gt.0) then
        mch=2
      else
        mch=1
      endif
      do ich=1,mch
        iu1=kch(ich)
        write(iu1,'('' RL max entropy, floating defaults'')')
        if (ih.eq.0) then
          txt='  (log likelihood)'
        else
          txt='  (rms likelihood)'
        endif
        write(iu1,'('' ih    '',i5,a)')ih,txt
        if (iw.eq.0)  then
          txt='  (no error bars read in)'
        else
          txt='  (error bars from input data)'
        endif
        write(iu1,'('' iw    '',i5,a)')iw,txt
        write(iu1,'('' pb0,pb1'',2f6.3)')pb0,pb1
        write(iu1,'('' ns    '',i5)')ns
        write(iu1,'('' ac    '',1pe10.2)')ac
        write(iu1,'('' nim   '',i5)')nim
        write(iu1,'('' al0, alf, nal   '',2f7.4,i4)')al0,alf,nal
        write(iu1,'('' clim  '',f7.4)')clim
        write(iu1,'('' ipri  '',2i3)')ipri
        write(iu1,'('' norm  '',2i3)')norm
        write(iu1,'('' wid, af '',e10.2,f7.4)')wid,af
      enddo
      pb0=pb0*pi2
      pb1=pb1*pi2
      return
      end

      subroutine fproj(pt,nil,nol,x,y)
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension x(*),y(*)
      real*4 pt(nri,nd)
      common /pdop/ indpt(nri,nd),indp(ndi,nr),ninpt(nd),ninp(nr)
      do i=1,nil
        y(i)=0
        m=ninpt(i)
        do j=1,m
          ict=indpt(j,i)
          y(i)=y(i)+pt(j,i)*x(ict)
        enddo
      enddo
      return
      end

      subroutine bproj(p,nil,nol,x,y)
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension x(*),y(*)
      real*4 p(ndi,nr)
      common /pdop/ indpt(nri,nd),indp(ndi,nr),ninpt(nd),ninp(nr)
      do j=1,nol
        y(j)=0
        n=ninp(j)
        do i=1,n
          ic=indp(i,j)
          y(j)=y(j)+p(i,j)*x(ic)
        enddo
      enddo
      return
      end

      subroutine prom(va,nv,pha,nph,vpa,nvp,pt,p)
      use dop_science
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension pha(*),vpa(*),vpl(nvpm+1),ww(nvpm)
      common/phases/dpha(npm)
      real*4 pt(nri,nd),p(ndi,nr)
      common /pdop/ indpt(nri,nd),indp(ndi,nr),ninpt(nd),ninp(nr)
      common /foldn/ nvf,indx(nr),indy(nr),nol
      common /cabsf/ wid,af
      real*8 kernel_cdf
      integer edge_bin
      va=max(abs(vpa(1)),abs(vpa(nvp)))
      dv=2*va/(nv-1)
      dvi=minval(vpa(2:nvp)-vpa(1:nvp-1))
      if (dvi.le.0.d0.or.nv.lt.2) stop 1
      do iv=1,nvp-1
        vpl(iv+1)=(vpa(iv)+vpa(iv+1))/2.d0
      enddo
      vpl(1)=vpa(1)-(vpa(2)-vpa(1))/2.d0
      vpl(nvp+1)=vpa(nvp)+(vpa(nvp)-vpa(nvp-1))/2.d0
      ninp=0
      ninpt=0
      ir=0
      icpm=0
      icm=0
      do ix=1,nv
        vx=dv*(ix-(nv+1)/2.d0)
        do iy=1,nv
          vy=dv*(iy-(nv+1)/2.d0)
          if (vx**2+vy**2.ge.(va+dv/10.d0)**2) cycle
          ir=ir+1
          if (ir.gt.nr) stop 1
          indx(ir)=ix
          indy(ir)=iy
          do ip=1,nph
            ww=0.d0
c Midpoint quadrature resolves the maximum velocity excursion to <=1/8 bin.
            nq=max(1,ceiling(abs(dpha(ip))*sqrt(vx*vx+vy*vy)/
     *                      min(dv,dvi)*8.d0))
            nq=max(nq,phase_subsamples)
            if (nq.gt.100000) then
              write(*,*) 'Exposure integration requires too many samples.'
              stop 1
            endif
            do iq=1,nq
              ph=pha(ip)+dpha(ip)*((iq-0.5d0)/nq-0.5d0)
              vr=vx*sin(ph)-vy*cos(ph)
              span=dv/2.d0+8.d0*instrument_sigma
              ki=max(1,edge_bin(vr-span,vpl,nvp+1))
              ka=min(nvp,edge_bin(vr+span,vpl,nvp+1))
              vfa=1.d0+af*exp(-min((vr/wid)**2,50.d0))
              if (vfa.le.0.d0) stop 1
              do iv=ki,ka
                w=kernel_cdf(vpl(iv+1)-vr,dv/2.d0,
     *                        instrument_sigma)-
     *            kernel_cdf(vpl(iv)-vr,dv/2.d0,
     *                        instrument_sigma)
                ww(iv)=ww(iv)+max(w,0.d0)/(nq*nph*vfa)
              enddo
            enddo
            do iv=1,nvp
              if (ww(iv).le.0.d0) cycle
              inp=ip+(iv-1)*nph
              icp=ninp(ir)+1
              ic=ninpt(inp)+1
              ninp(ir)=icp
              ninpt(inp)=ic
              icpm=max(icpm,icp)
              icm=max(icm,ic)
              if (icp.le.ndi.and.ic.le.nri) then
                p(icp,ir)=ww(iv)
                indp(icp,ir)=inp
                pt(ic,inp)=ww(iv)
                indpt(ic,inp)=ir
              endif
            enddo
          enddo
        enddo
      enddo
      nol=ir
      if (icpm.gt.ndi.or.icm.gt.nri) then
        write(10,*) 'ndi,nri:',ndi,nri
        write(10,*) 'projection matrix array too small. Need ndi,nri=',
     *              icpm,icm
        call wpar(npm,nvpm,nvm,icm,icpm)
        stop 2
      endif
      end

      integer function edge_bin(x,edges,n)
      implicit none
      integer n,lo,hi,mid
      real*8 x,edges(n)
      lo=0
      hi=n
      do while (hi-lo.gt.1)
        mid=(lo+hi)/2
        if (edges(mid).le.x) then
          lo=mid
        else
          hi=mid
        endif
      enddo
      edge_bin=lo
      end

      real*8 function normal_second(x,sigma)
      implicit none
      real*8 x,sigma,z,cdf,pdf
      z=x/sigma
      cdf=0.5d0*erfc(-z/sqrt(2.d0))
      pdf=exp(-0.5d0*z*z)/sqrt(8.d0*atan(1.d0))
      normal_second=0.5d0*((x*x+sigma*sigma)*cdf+x*sigma*pdf)
      end

      real*8 function kernel_cdf(x,a,sigma)
c CDF of the normalized triangular pixel basis convolved with a Gaussian LSF.
      implicit none
      real*8 x,a,sigma,dist,normal_second,z
      if (sigma.le.0.d0) then
        kernel_cdf=dist((x+a)/(2.d0*a))
      else if (x.le.-a-8.d0*sigma) then
        kernel_cdf=0.d0
      else if (x.ge.a+8.d0*sigma) then
        kernel_cdf=1.d0
      else
        z=(normal_second(x+a,sigma)-2.d0*normal_second(x,sigma)+
     *     normal_second(x-a,sigma))/(a*a)
        kernel_cdf=max(0.d0,min(1.d0,z))
      endif
      end

      real*8 function dist(x)
      implicit real*8 (a-h,o-z)
      if (x.lt.0.) then
        dist=0.
      else
        if (x.lt.0.5) then
          dist=2*x**2
        else
          if (x.lt.1.) then
            dist=1-2*(1.-x)**2
          else
            dist=1.
          endif
        endif
      endif
      return
      end

      subroutine erri(s,n1,n2)
      implicit real*8 (a-h,o-z)
      character*(*) s
      dimension x(2)
      write(*,'(a,2i10)')s,n1,n2
      write(10,'(a,2i10)')s,n1,n2
      stop 1
c      return
      end

      subroutine errb(s,x1,x2)
      implicit real*8 (a-h,o-z)
      character*(*) s
      dimension x(2)
      write(*,'(a,1p2e12.4)')s,x1,x2
      stop 1
c      return
      end

      subroutine weighd(f,f0,pha,nph,nvp,eb,iw)
      use dop_science
      implicit real*8 (a-h,o-z)
      dimension f(*),f0(*),pha(*),eb(*)
      do ip=1,nph
        do iv=1,nvp
          k=ip+(iv-1)*nph
          stat_weights(k)=0.d0
          if (.not.phase_valid(ip)) cycle
          stat_weights(k)=1.d0
          if (iw.eq.1) stat_weights(k)=1.d0/eb(k)**2
        enddo
      enddo
      end

      subroutine rldop(pht,p,pt,al,ac,flim,nim,ps,nil,nol,nit,hs,
     *           ns,ih,ipri,s,ierr)
      implicit real*8 (a-h,o-z)
      dimension pht(*),ps(*),epc(1)
      integer inc(1,1),nsm(1,1),ico(1)
      real*4 p(*),pt(*)
      inc=1
      nsm=ns
      ico=1
      epc=1.d0
      call rlfldn(pht,p,pt,al,ac,flim,nim,ps,nil,nol,nit,hs,
     *           1,inc,nsm,epc,ico,0.d0,1,ih,ipri,0.d0,s,ierr)
      end

      subroutine chims(ps,chi,ns)
c default map routine to interface dopmap with emap routines
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension ps(*),chi(nr)
      dimension x2(nvm,nvm),y2(nvm,nvm)
      common /foldn/ nv,indx(nr),indy(nr),nol
c reset x2
      do i=1,nv
        do j=1,nv
          x2(i,j)=0
        enddo
      enddo
c fold x into x2 using `backdoor' information from /foldn/
      do j=1,nol
        x2(indx(j),indy(j))=ps(j)
      enddo
      ic=0
      jc=0
      im=1
      call chimg(x2,y2,nv,ic,nv,jc,im,ns)
c unfold y2
      do j=1,nol
        chi(j)=y2(indx(j),indy(j))
      enddo
      return
      end

      subroutine wpar(npm,nvpm,nvm,nri,ndi)
      implicit none
      integer npm,nvpm,nvm,nri,ndi
      open(11,file='emap.par.new',status='replace')
      write(11,*) '      parameter (npm=',npm,',nvpm=',nvpm,
     *            ',nvm=',nvm,')'
      write(11,*) '      parameter (nd=npm*nvpm,nr=nvm*nvm)'
      write(11,*) '      parameter (nri=',nri,',ndi=',ndi,')'
      write(11,*) '      parameter (nf=nd,nfp=ndi,nrp=nri)'
      write(11,*) '      parameter (ni=nvm,nj=nvm)'
      write(11,*) '      parameter (nch=1,nsu=1)'
      close(11)
      end

      subroutine chimg(x2,y2,nam,ic,nbm,jc,im,ns)
c gaussian smearing with half-width ns, by consecutive 1-d smearing in x and y
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      parameter (ms=100)
      dimension x2(ni,nj), y2(ni,nj)
      real*8 smp(-ms:ms)
      integer, save :: ns0=-1,ns1=0
      save smp
      dimension t1(-ms:ni+ms),t2(-ms:nj+ms)
      if (ns0.ne.ns) then
c (re)make smearing matrix smp if value of ns is new.
        w=1.5
c actual smearing width used is w times 1/e width
        ns1=min(int(w*ns),ms,nam-1,nbm-1)
        if (ns1.gt.ms) call erri('ns1>ms in chimg',ns1,ms)
        ns0=ns
        s=0
        do i=-ns1,ns1
          smp(i)=exp(-float(i*i)/ns**2)
          s=s+smp(i)
        enddo
c normalize
        do i=-ns1,ns1
          smp(i)=smp(i)/s
        enddo
      endif

c smear in j-direction
      do i=1,nam
c extend ith column of x2 into t2, taking properties of boundary into account
        do j=-ns1,0
          if (jc.eq.1) then
c wraparound
            t2(j)=x2(i,j+nbm)
          else
c fold back
            t2(j)=x2(i,2-j-1)
          endif
        enddo
        do j=1,nbm
          t2(j)=x2(i,j)
        enddo
        do j=nbm+1,nbm+ns1
          if (jc.eq.1) then
             t2(j)=x2(i,j-nbm)
          else
             t2(j)=x2(i,2*nbm-j+1)
          endif
        enddo
c smear t2, result into y2
        do j=1,nbm
          s=0
          do js=-ns1,ns1
             s=s+t2(j+js)*smp(js)
          enddo
          y2(i,j)=s
        enddo
      enddo

c smear in i-direction
      do j=1,nbm
c extend jth row of y2 into t1, taking properties of boundary into account
        do i=-ns1,0
          if (ic.eq.1) then
c wraparound
            t1(i)=y2(i+nam,j)
          else
c fold back
            t1(i)=y2(2-i-1,j)
          endif
        enddo
        do i=1,nam
          t1(i)=y2(i,j)
        enddo
        do i=nam+1,nam+ns1
          if (ic.eq.1) then
             t1(i)=y2(i-nam,j)
          else
             t1(i)=y2(2*nam-i+1,j)
          endif
        enddo
c smear t1, result into y2
        do i=1,nam
          s=0
          do is=-ns1,ns1
             s=s+t1(i+is)*smp(is)
          enddo
          y2(i,j)=s
        enddo
      enddo
      return
      end

      subroutine rlfldn(pht,p,pt,al,ac,flim,nim,ps,nil,nol,nit,hs,
     *           nc,inc,ns,epc,ico,rdf,msu,ih,ipri,wbkg,s,ierr)
c Positive, preconditioned gradient ascent with Newton acceleration and
c an Armijo line search. No imposed flux sum: masked data have no effect.
      use dop_science, only: solver_objective,solver_stationarity
      use, intrinsic :: ieee_arithmetic
      implicit real*8 (a-h,o-z)
      real*4 dpx
      include 'emap.par'
      dimension pht(*),ps(*),inc(nch,nsu),ns(nch,nsu),epc(*),ico(*)
      real*4 p(nfp,nr),pt(nrp,nf)
      dimension ph(nf),dha(nr),dsa(nr),chi(nr,nch),tem(nr),ex(nr)
      dimension dps(nr),dps0(nr),up(nr),trial(nr),grad(nr)
      common /difpix/ dpx(nr,4)
      ierr=1
      dps0=0.d0
      dpx=0.
      s=0.d0
      hs=0.d0
      do nit=0,nim
        call fproj(pt,nil,nol,ps,ph)
        call like(pht,ph,p,nil,nol,ih,h,dha)
        dsa=0.d0
        if (al.gt.0.d0) call entr(ps,nil,nol,s,dsa,chi,tem,ex,
     *      nc,inc,ns,epc,ico,rdf,msu,wbkg)
        hs=h+al*s
        solver_objective=hs
        if (.not.ieee_is_finite(hs)) then
          ierr=3
          return
        endif
        su=0.d0
        sh=0.d0
        ss=0.d0
        do j=1,nol
          grad(j)=dha(j)+al*dsa(j)
          dps(j)=ps(j)*grad(j)
          su=su+dps(j)**2
          sh=sh+(ps(j)*dha(j))**2
          ss=ss+(ps(j)*al*dsa(j))**2
        enddo
        denom=max(sqrt(sh)+sqrt(ss),1.d-12)
        del=sqrt(su)/denom
        solver_stationarity=del
        if (nit.eq.0.or.mod(nit,10).eq.0.or.del.le.ac)
     *      write(10,*) 'iteration objective stationarity ',nit,hs,del
        if (del.le.ac.or.sqrt(su).le.1.d-12) then
          ierr=0
          return
        endif
        if (nit.eq.nim) return
        if (nit.eq.0) dps0(1:nol)=dps(1:nol)
        call lamun(pht,ps,dps0,dps,pt,ph,chi,tem,ex,al,nil,nol,
     *             fla,fmu,nc,inc,ns,epc,ico,rdf,msu,ih,wbkg)
        up(1:nol)=fla*dps0(1:nol)+fmu*dps(1:nol)
        slope=dot_product(grad(1:nol),up(1:nol))
        if (.not.ieee_is_finite(slope).or.slope.le.0.d0) then
          up(1:nol)=dps(1:nol)
          slope=dot_product(grad(1:nol),up(1:nol))
        endif
        step=1.d0
        do j=1,nol
          if (up(j).lt.0.d0) step=min(step,-0.99d0*ps(j)/up(j))
        enddo
        do ils=1,60
          trial(1:nol)=ps(1:nol)+step*up(1:nol)
          call fproj(pt,nil,nol,trial,ph)
          call like(pht,ph,p,nil,nol,ih,hnew,dha)
          snew=0.d0
          if (al.gt.0.d0) call entr(trial,nil,nol,snew,dsa,
     *        chi,tem,ex,nc,inc,ns,epc,ico,rdf,msu,wbkg)
          qnew=hnew+al*snew
          if (ieee_is_finite(qnew).and.
     *        qnew.ge.hs+1.d-4*step*slope) exit
          step=step/2.d0
        enddo
        if (ils.gt.60) then
          ierr=2
          return
        endif
        ps(1:nol)=trial(1:nol)
        dps0(1:nol)=step*up(1:nol)
        do j=1,nol
          dpx(j,mod(nit,4)+1)=grad(j)
        enddo
      enddo
      end

      subroutine like(pht,ph,p,nil,nol,ih,h,dha)
      use dop_science
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension pht(*),ph(*),dha(*),res(nf),phh(nf)
      real*4 p(nfp,nr)
      h=0.d0
      if (ih.eq.0.and.stat_iw.eq.0) then
        do i=1,nil
          phh(i)=0.d0
          if (stat_weights(i).eq.0.d0) cycle
          pp=max(ph(i),tiny(1.d0))
          h=h+pht(i)*log(pp)-pp
          phh(i)=pht(i)/pp-1.d0
        enddo
      else
        do i=1,nil
          res(i)=pht(i)-ph(i)
        enddo
        call apply_precision(res,phh,q)
        h=-0.5d0*q
      endif
      call bproj(p,nil,nol,phh,dha)
      end

      subroutine entr(ps,nil,nol,s,sp,chi,tem,ex,nc,inc,ns,epc,ico,
     *                rdf,msu,wbkg)
c entropy part of `objective function' Q.
c input: ps (current reconstruction),
c output: s (S), sp (dS/dpsi)
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension ps(*),sp(*),chi(nr,nch),tem(*),ex(*),inc(nch,nsu),
     *          ns(nch,nsu),epc(*),ico(*)
      dimension x(nr),y(nr,nch)
      if (nc.gt.nch) call erri('nc>nch in entr',nc,nch)
c make chi's
      call chims(ps,chi(:,1),ns(1,1))
c entropy
      s=0
      do j=1,nol
        tem(j)=log(ps(j))
        ex(j)=0
        do k=1,nc
          if (epc(k).ne.0.) then
            fe=epc(k)*log(chi(j,k))
            ex(j)=ex(j)+fe
            tem(j)=tem(j)-fe
          endif
        enddo
        ex(j)=exp(ex(j))
        s=s-ps(j)*(tem(j)-1)-ex(j)
      enddo
c derivative wr to psi_j, first the psilogpsi part
      do j=1,nol
        sp(j)=-tem(j)
      enddo
c add terms due to dependence of chi on psi
      do k=1,nc
        if (epc(k).ne.0.) then
          do j=1,nol
            x(j)=(ps(j)-ex(j))/chi(j,k)
          enddo
          call chims(x,y(:,k),ns(k,1))
          do j=1,nol
            sp(j)=sp(j)+epc(k)*y(j,k)
          enddo
        endif
      enddo
      return
      end

      subroutine lamun(pht,ps,dps1,dps2,pt,ph,chi,tem,ex,al,nil,nol,
     *                 fla,fmu,nc,inc,ns,epc,ico,rdf,msu,ih,wbkg)
c optimal lambda and mu from first and second derivatives of Q
      use dop_science
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension pht(*),ps(*),dps1(*),dps2(*),ph(*),chi(nr,nch),
     *          tem(*),ex(*),inc(nch,nsu),ns(nch,nsu),epc(*),ico(*)
      real*4 pt(nrp,nf)
      dimension dph1(nf),dph2(nf),dch1(nr,nch),dch2(nr,nch)
      dimension rw(nf),z1(nf),z2(nf),res(nf)
c changes in phi and chi corresponding to dps0, dps
      call fproj(pt,nil,nol,dps1,dph1)
      call fproj(pt,nil,nol,dps2,dph2)
      if (al.gt.0.) then
        call chims(dps1,dch1(:,1),ns(1,1))
        call chims(dps2,dch2(:,1),ns(1,1))
      endif
c first and second derivatives
      hl=0
      hll=0
      hm=0
      hlm=0
      hmm=0
      sl=0
      sll=0
      sm=0
      slm=0
      smm=0
c Exact first and second directional derivatives of the likelihood.
      if (ih.eq.0.and.stat_iw.eq.0) then
        do i=1,nil
          if (stat_weights(i).eq.0.d0) cycle
          pp=max(ph(i),tiny(1.d0))
          dpp=pht(i)/pp-1.d0
          xl=dph1(i)/pp
          xm=dph2(i)/pp
          hl=hl+dph1(i)*dpp
          hll=hll-xl**2*pht(i)
          hlm=hlm-xl*xm*pht(i)
          hm=hm+dph2(i)*dpp
          hmm=hmm-xm**2*pht(i)
        enddo
      else
        do i=1,nil
          res(i)=pht(i)-ph(i)
        enddo
        call apply_precision(res,rw,q)
        call apply_precision(dph1,z1,q1)
        call apply_precision(dph2,z2,q2)
        do i=1,nil
          hl=hl+dph1(i)*rw(i)
          hm=hm+dph2(i)*rw(i)
          hlm=hlm-dph1(i)*z2(i)
        enddo
        hll=-q1
        hmm=-q2
      endif
c entropy terms
      if (al.gt.0.) then
      do j=1,nol
        pl=dps1(j)/ps(j)
        pm=dps2(j)/ps(j)
        tl=0
        tm=0
        tll=0
        tmm=0
        tlm=0
        do k=1,nc
          if (epc(k).ne.0.) then
            dc1=dch1(j,k)/chi(j,k)
            dc2=dch2(j,k)/chi(j,k)
            tl=tl+epc(k)*dc1
            tm=tm+epc(k)*dc2
            tll=tll+epc(k)*dc1**2
            tmm=tmm+epc(k)*dc2**2
            tlm=tlm+epc(k)*dc1*dc2
          endif
        enddo
        pz=ps(j)-ex(j)
        sl=sl-dps1(j)*tem(j)+pz*tl
        sm=sm-dps2(j)*tem(j)+pz*tm
        sll=sll-ps(j)*(pl-tl)**2+pz*(tl**2-tll)
        slm=slm-ps(j)*(pl-tl)*(pm-tm)+pz*(tl*tm-tlm)
        smm=smm-ps(j)*(pm-tm)**2+pz*(tm**2-tmm)
      enddo
      endif

      ql=hl+al*sl
      qll=hll+al*sll
      qlm=hlm+al*slm
      qm=hm+al*sm
      qmm=hmm+al*smm
c Robust two-direction Newton proposal; the caller applies line search.
      dd=qmm*qll-qlm**2
      denom=qll**2+qmm**2
      fla=0.d0
      fmu=1.d0
      if (denom.le.tiny(1.d0)) return
      dr=abs(dd)/denom
      if (dr.lt.1.d-8.or.dd.le.0.d0.or.qll.ge.0.d0) then
        if (qmm.lt.-tiny(1.d0)) fmu=-qm/qmm
      else
        fmu=(qlm*ql-qll*qm)/dd
        fla=(qlm*qm-qmm*ql)/dd
      endif
      return
      end

      function rrms(f1,eb,iw,f0,n)
      use dop_science
      implicit real*8 (a-h,o-z)
      dimension f1(*),eb(*),f0(*),res(n),weighted(n)
      do i=1,n
        res(i)=f1(i)-f0(i)
      enddo
      call apply_precision(res,weighted,q)
      rrms=sqrt(max(q,0.d0)/statistic_rank)
      end
