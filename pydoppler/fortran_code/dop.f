      program dopmap
c Version 2.0
c HCS 3-3-1999
c version for arbitrary phase bins
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
      pi2=8*atan(1.)

c read parameter file
      call repar(ns,ac,nim,
     *     al0,alf,nal,clim,ipri,spec,ih,iw,pb0,pb1,norm)
c set baseline level to be added to spectrum (fraction of the avg of entire spectrum)
      bfa=0.1
      write(*,'('' HOLAQQ NOW 2Q'')')
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
      if (iw.eq.1) read(4,*) (eb(i),i=1,nil)
      close(4)
      call clock(t1)
      call prom(va,nv,pha,nph,vp,nvp,pt,p)
      call clock(t2)
      write(10,'('' cpu for geometry'',f10.2)')t2-t1
      if (ipri.gt.0) write(*,'('' cpu for geometry'',f10.2)')t2-t1
      t1=t2
c compute a baseline flux to add to dat
      call basel(dat,nph,nvp,bfa,bf,nv,va,pt,datb,vmb)
c add baseline to dat, put weighting map into dato
      do i=1,nil
        dat(i)=dat(i)+datb(i)
        if (dat(i).lt.bf/3) then
          dat(i)=1e-5*bf
          dato(i)=0
        else
          dato(i)=1
        endif
      enddo
c normalize to constant wavelength-integrated flux
      if (norm.eq.1) call norms(dat,pb0,pb1,pha,nph,nvp)
c sum of light
      sdat=0
      do i=1,nil
        sdat=sdat+dat(i)
      enddo

c make error bars relative
      do i=1,nil
        eb(i)=eb(i)/dat(i)
      enddo
c weights
      call weighd(dat,datb,pha,nph,nvp,eb,iw)

c reconstruction
c initialize (flat)
      write(*,'('' HOLAQQQQ'')')
      do j=1,nol
        vmo(j)=sdat/nol
      enddo
c decreasing sequence of alfas until required fit to input data
      call alfs(dat,eb,iw,p,pt,al0,alf,nal,ac,nim,vmo,
     *           dato,nil,nol,nit,al,rr,ns,ih,clim,ipri,ier,t2,sx)
      if (ier.eq.2) then
        write(*,'('' no converged models, returning starting image'')')
        al=al0
        rr=100
      endif
      if (ier.eq.1) write(*,'('' last converged model has rr='',
     *    f7.4,'', >clim='',f7.4)') rr,clim
      if (ipri.gt.0) write(*,'('' entropy'', 1pe12.4)')sx
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
      fident=sin(100*(s1/s2))
c normalize dopmap to input units (probably mJy's) x Hz per (cm/s)^2 in v-space
c velocity grid size
      dv=2*va/(nv-1)
      delv2=dv**2
c frequency grid size
      dev=(vp(nvp)-vp(1))/(nvp-1)
      denu=dev/w0*1e8
      do i=1,nol
        vmo(i)=(vmo(i)-vmb(i))/delv2*denu
      enddo
c make unfolded reconstructed image
      do j=1,nol
        reco(indx(j),indy(j))=vmo(j)
      enddo

      open(3,file='dop.out')
c write array sizes, w0, identification
      write(3,'(3i5,f10.3,f20.15)')nph,nvp,nvm,w0,fident
c copy second line of input data file
      write(3,*)idat
c copy phase and velocity coordinates of input spectrum
      write(3,'(1p8e13.4E3)')(pha(i),i=1,nph)
      write(3,*)iph
      if (iph.eq.1) write(3,'(1p8e13.4E3)')(dpha(i),i=1,nph)
      write(3,'(1p8e13.4E3)')(vp(i),i=1,nvp)
c write input spectrum
      write(3,'(1p8e13.4E3)')(dat(i)-datb(i),i=1,nil)
      write(3,'(2i3,2f7.3,i4,2e13.4E3,f7.3,i3,e13.4E3,f7.3)')
     *      ih,iw,pb0/pi2,pb1/pi2,ns,ac,al,clim,norm,wid,af
c write dopmap and reco
      write(3,'(i4,1pe13.4E3,'' dopmap '')')nv,va
      write(3,'(1p6e13.4E3)')((reco(i,j),j=1,nv),i=1,nv)
      write(3,'(i6,'' reconstructed data '')')nil
      write(3,'(1p6e13.4E3)')(dato(i)-datb(i),i=1,nil)
      write(3,'(i6,'' 4  change maps '')')nv
      do ii=1,4
        do j=1,nol
          dpxu(indx(j),indy(j))=dpx(j,ii)
        enddo
        write(3,'(1p6e13.4E3)')((dpxu(i,j),j=1,nv),i=1,nv)
      enddo
      close(3)
      end

c      block data
c      include '/usr/include/fpdc.h'
c      include '/usr/include/fpdt.h'
c      end

      subroutine basel(dat,nph,nvp,bfa,bf,nv,va,pt,datb,vmb)
c compute a baseline flux to get a flat background level in dopmap
c datb is the data corresponding to a flat (circular) dopmap
c vmb is the corresponding flat dopmap
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      common /foldn/ nvf,indx(nr),indy(nr),nol
      dimension dat(*),datb(*),vmb(*)
      if (nol.gt.nr) call errb('nol>nr',nol,nr)

c make circle of radius va on a zero background
      fv2=nv/2.
      dv=2*va/(nv-1)
      do ir=1,nol
        ix=indx(ir)
        iy=indy(ir)
        vx=dv*(ix-fv2-0.5)
        vy=dv*(iy-fv2-0.5)
        v=sqrt(vx**2+vy**2)
        vmb(ir)=0.
        if (v.lt.va) vmb(ir)=1
      enddo

      nil=nph*nvp
c produce spectrum from this map
      call fproj(pt,nil,nol,vmb,datb)
c value at v=0
      d0=datb(nvp/2*nph)
c find avg of dat
      avg=0.
      do i=1,nil
        avg=avg+dat(i)
      enddo
      avg=avg/nil
c baseline level to add
      bf=bfa*avg
c normalize datb
      do i=1,nil
        datb(i)=bf*datb(i)/d0
      enddo
c normalize background map
      do j=1,nol
        vmb(j)=vmb(j)*bf/d0
      enddo
      return
      end

      subroutine blanks(pha,nph,pb0,pb1)
c make blanking array corresponding to (pb0,pb1)
c pb0-pb1: phase range with weight 0 (used for ignoring eclipsed part of data)
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension pha(*)
      common/blank/bl(npm)
      logical bl
      pi2=8*atan(1.)
c get max and min values of phase
      phmin=pha(1)
      phmax=pha(1)
      do i=2,nph
        if (pha(i).lt.phmin) phmin=pha(i)
        if (pha(i).gt.phmax) phmax=pha(i)
      enddo
c number of cycles present
      nphmi=int(phmin/pi2)
      if (nphmi.lt.0.) nphmi=nphmi-1
      nphma=int(phmax/pi2)
      if (nphma.lt.0.) nphma=nphma-1
c reduce blanking region to phases around 0
      npb1=int(pb1/pi2)
      if (npb1.lt.0.) npb1=npb1-1
      pb0=pb0-npb1*pi2
      pb1=pb1-npb1*pi2
      do i=1,nph
c use only non-blanked region for making average
        bl(i)=.false.
        do ic=nphmi,nphma+1
          if (pha(i).gt.pb0+ic*pi2.and.pha(i).lt.pb1+ic*pi2)
     *      bl(i)=.true.
        enddo
      enddo
      return
      end

      subroutine norms(dat,pb0,pb1,pha,nph,nvp)
c normalize by a phase-dependent factor to make wavelength-integrated
c light constant with phase.
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension dat(npm,nvpm),pha(*),wi(npm)
      common/blank/bl(npm)
      logical bl
      pi2=8*atan(1.)
      ws=0
      np=0
c use only non-blanked region for making average
      do i=1,nph
        wi(i)=0
        if (.not.bl(i)) then
          np=np+1
          do j=1,nvp
            wi(i)=wi(i)+dat(i,j)
          enddo
        endif
        ws=ws+wi(i)
      enddo
      ws=ws/np
      do i=1,nph
        wi(i)=wi(i)/ws
      enddo
c measure minimum of wi
      wim=0
      do i=1,nph
        if (wi(i).lt.wim) wim=wi(i)
      enddo
c set an overall level to be added before normalization
      wl=0
      wlm=0.3
      if (wim.lt.wlm) then wl=wlm-wim
      do i=1,nph
        if (wi(i).ne.0.) then
          do j=1,nvp
            dat(i,j)=dat(i,j)*(1+wl)/(wi(i)+wl)
          enddo
        endif
      enddo
      return
      end

      subroutine alfs(dat,eb,iw,p,pt,al0,alf,nal,ac,nim,vmo,
     *           dato,nil,nol,nit,alx,rr,ns,ih,clim,ipri,ier,t2,sx)
c sequence of models with decreasing alfas, starting at al0, in factors
c of alf, max no of models nal, each converged to accuracy ac, until
c rms difference with input data dat is less than clim
c ier=0:  normal return
c ier=1:  last converged model differs more than clim from input light curve
c ier=2:  no converged models in entire sequence
      implicit real*8 (a-h,o-z)
      include 'emap.par'
c dopmaps
      dimension vm1(nr),vmo(*)
c projection matrix
      real*4 p(nfp,nr),pt(nrp,nf)
c input data
      dimension dat(nd),eb(nd),dat1(nd),dato(nd)
c parameters of penalty functions
      logical cont
      ier=0
c copy input into temporary image
      do j=1,nol
        vm1(j)=vmo(j)
      enddo
c initialize number of converged models
      nco=0
c alpha- loop
      rr=100
      ial=0
      cont=.true.
      call clock(t2)
      t1=t2
c if nal<0 converge abs(nal) models irrespective of clim
      if (nal.lt.0) then
        nal1=-nal
      else
        nal1=nal
      endif
      do while (ial.lt.nal1.and.(rr.gt.clim.or.nal.lt.0).and.cont)
        ial=ial+1
c next alfa value
        if (nal.lt.0) then
          al=al0*alf**(1-ial)
        else
          al=alfa(al0,ial,alf,rr,clim)
        endif
        aci=ac
        flim=1e5
c converge with Lucy's ME scheme
        call rldop(dat,p,pt,al,ac,flim,nim,vm1,nil,nol,nit,hs,
     *           ns,ih,ipri,sf,ierr)
        if (ierr.eq.1) write(10,'('' not converged'')')
        if (ierr.eq.2) write(10,'('' not monotonic'')')
        if (ierr.eq.3) write(10,'('' not converged, not monotonic'')')
        if (ierr.eq.4) write(10,'('' stuck at vanishing updates'')')
        if (ipri.gt.0) then
          if (ierr.eq.1) write(*,'('' not converged'')')
          if (ierr.eq.2) write(*,'('' not monotonic'')')
          if (ierr.eq.3) write(*,'('' not converged, not monotonic'')')
          if (ierr.eq.4) write(*,'('' stuck at vanishing updates'')')
        endif
c reconstructed spectrum
        call fproj(pt,nil,nol,vm1,dat1)
c converged at this alpha, save.
          nco=nco+1
          do j=1,nol
            vmo(j)=vm1(j)
          enddo
          do k=1,nil
            dato(k)=dat1(k)
          enddo
          alx=al
          sx=sf
c rms difference w/r input light curve
        rr0=rr
        rr=rrms(dat,eb,iw,dat1,nil)
        if (iw.eq.0) then
c measure accuracy of reconstruction relative to point-to-point noise in data
          rn=rmsn(dat,dat1)
          rr=rr/rn
        endif
c summary of this iteration to activity log file
        write(10,'('' ni, al, hs, rr:'',i5,
     *       f8.5,e17.9,2f8.5)')nit,al,hs,rr,clim
        call clock(t2)
        write(10,'('' cpu for iteration'',f8.2)')t2-t1
        if (ipri.gt.0) then
          write(*,'('' ni, al, hs, rr:'',i5,
     *       f8.5,e17.9,2f8.5)')nit,al,hs,rr,clim
          write(*,'('' cpu for iteration'',f8.2)')t2-t1
        endif
        t1=t2
c continue if no converged models found yet, or if acceptable error
c return from last model
        cont = nco.eq.0.or.ierr.eq.0.or.ierr.eq.2
      enddo
      if (nco.eq.0) ier=2
      if (.not.cont) then
        ier=1
        rr=rr0
      endif
      return
      end

      subroutine repar(ns,ac,nim,
     *     al0,alf,nal,clim,ipri,spec,ih,iw,pb0,pb1,norm)
c read input parameters
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      character*(*) spec
      character*40 txt
c output channels
      dimension kch(2)
      common /fold/ nvx,nvy
      common /cabsf/ wid,af
      pi2=8*atan(1.)
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
c projection matrix
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension pha(*),vpa(*)
      common/phases/dpha(npm)
c pha: mid-exposure phases, dpha: length of phase intervals
      real*4 pt(nri,nd),p(ndi,nr)
c ninpt: stores length of first index of pt, ninp: same for p
      common /pdop/ indpt(nri,nd),indp(ndi,nr),ninpt(nd),ninp(nr)
      common /foldn/ nvf,indx(nr),indy(nr),nol
c cos/sin tables
      dimension cfl(npm),sfl(npm),cfu(npm),sfu(npm)
c lower and upper edges phase bins
      dimension phl(npm),phu(npm),dph(npm)
      dimension vpl(-1:nvpm+2)
c array for temp storage of sums of projection matrix
      dimension s(nr)
c parameters for `central absorption'
      common /cabsf/ wid,af
      va=dmax1(abs(vpa(1)),vpa(nvp))

      nil=nph*nvp
      fv2=nv/2.
      dv=2*va/(nv-1)
      dvi=(vpa(nvp)-vpa(1))/(nvp-1)
      do ip=1,nph
        phl(ip)=pha(ip)-dpha(ip)/2
        phu(ip)=pha(ip)+dpha(ip)/2
        cfl(ip)=cos(phl(ip))
        cfu(ip)=cos(phu(ip))
        sfl(ip)=sin(phl(ip))
        sfu(ip)=sin(phu(ip))
      enddo
c lower edges of velocity bins
      do iv=1,nvp-1
        vpl(iv+1)=(vpa(iv)+vpa(iv+1))/2
      enddo
      vpl(1)=vpa(1)-(vpa(2)-vpa(1))/2
      vpl(nvp+1)=vpa(nvp)+(vpa(nvp)-vpa(nvp-1))/2
      vpl(0)=vpl(1)-(vpl(2)-vpl(1))
      vpl(-1)=-10*va
      vpl(nvp+2)=10*va
c reset index arrays
      do i=1,nd
        ninpt(i)=0
      enddo
      do i=1,nr
        ninp(i)=0
      enddo
c reset x-y counter and bounds check flag
      ir=0
      ierp=0
      do ix=1,nv
c x-velocity (grid defined by this statement). Min and max of vx are (-va,va)
        vx=dv*(ix-fv2-0.5)
        do iy=1,nv
c y-velocity (grid defined by this statement). Min and max of vy are (-va,va)
          vy=dv*(iy-fv2-0.5)
c use only points inside circle of radius va
          if (vx**2+vy**2.lt.(va+dv/10)**2) then
            ir=ir+1
            if (ir.gt.nr) call erri('ir>nr in prom',ir,nr)
            indx(ir)=ix
            indy(ir)=iy
            do ip=1,nph
c line-of sight velocity at phase bin edges
              vp0=vx*sfl(ip)-vy*cfl(ip)
              vp1=vx*sfu(ip)-vy*cfu(ip)
              if (vp0.gt.vp1) then
                vpma=vp0+dv/2
                vpmi=vp1-dv/2
              else
                vpma=vp1+dv/2
                vpmi=vp0-dv/2
              endif
c estimate near which line-sight bins these are
              ka=(vpma-vpl(1))/dvi+1
              ki=(vpmi-vpl(1))/dvi+1
              if (vpmi.lt.vpl(-1))
     *            call errb('vpmi,vpl(-1):',vpmi,vpl(-1))
              if (vpma.gt.vpl(nvp+2))
     *            call errb('vpma,vpl(nvp+2):',vpma,vpl(nvp+2))
c get first vpl(i)<vpma
              ka=min(ka,nvpm+1)
              ki=max(0,ki)
              ki=min(ki,nvpm+1)
              ka=max(0,ka)
              if (vpl(ka).lt.vpma) then
                do while (vpl(ka+1).lt.vpma)
                  ka=ka+1
                enddo
              else
                do while (vpl(ka-1).gt.vpma)
                  ka=ka-1
                enddo
              endif
c get first vpl(i)<vpmi
              if (vpl(ki).lt.vpmi) then
                do while (vpl(ki+1).lt.vpmi)
                  ki=ki+1
                enddo
              else
                do while (vpl(ki-1).gt.vpmi)
                  ki=ki-1
                enddo
              endif
              ka=min(ka,nvpm)
              ki=max(1,ki)
              do iv=ki,ka
c compute index of this phase-velocity point in data array
                inp=ip+(iv-1)*nph
                x0=(vpl(iv)-vpmi)/(vpma-vpmi)
                x1=(vpl(iv+1)-vpmi)/(vpma-vpmi)
                w=dist(x1)-dist(x0)
c central absorption fudge
                vpar=(vpma+vpmi)/2
                vfa=1+af*exp(-dmin1((vpar/wid)**2,50.d0))
                w=w/vfa
c up counter for p
                icp=ninp(ir)+1
                if (icp.gt.ndi) then
                  ierp=2
                  icpm=max(icp,icpm)
                endif
                ninp(ir)=icp
                if (ierp.ne.2) then
c store matrix element and save data index
                  p(icp,ir)=w
                indp(icp,ir)=inp
                endif
c up counter for pt of number of hits in this bin
                ic=ninpt(inp)+1
                if (ic.gt.nri) then
                  ierp=1
                  icm=max(ic,icm)
                endif
                ninpt(inp)=ic
                if (ierp.ne.1) then
c store transpose element and save reconstruction index
                  pt(ic,inp)=w
                  indpt(ic,inp)=ir
                endif
              enddo
            enddo
          endif
        enddo
      enddo
      if (ierp.ne.0) then
        if (icpm.lt.ndi)icpm=ndi
        if (icm.lt.nri)icm=nri
        write(*,'('' ndi,nri:'',2i6)')ndi,nri
        write(*,'('' projection matrix array too small. Need ndi,nri='',
     *        2i6,2f7.3)')icpm,icm,float(icpm)/ndi,float(icm)/nri
        write(10,'('' ndi,nri:'',2i6)')ndi,nri
        write(10,'(''projection matrix array too small. Need ndi,nri='',
     *        2i6)')icpm,icm
        call wpar(npm,nvpm,nvm,icm,icpm)
        stop
      endif
      nol=ir
c normalize p,pt
      do ir=1,nol
        s(ir)=0
        do j=1,ninp(ir)
          s(ir)=s(ir)+p(j,ir)
        enddo
        do j=1,ninp(ir)
          p(j,ir)=p(j,ir)/s(ir)
        enddo
      enddo
      do inp=1,nil
        do j=1,ninpt(inp)
          pt(j,inp)=pt(j,inp)/s(indpt(j,inp))
        enddo
      enddo
c count storage requirements
      ima=0
      mninp=0
      do ir=1,nol
        len=ninp(ir)
        mninp=mninp+len
        if (len.gt.ima) ima=len
      enddo
      jma=0
      mninpt=0
      do inp=1,nil
        len=ninpt(inp)
        mninpt=mninpt+len
        if (len.gt.jma) jma=len
      enddo
c      write(*,'('' max of data, reco counters: '',2i7)')ima,jma
c      write(*,'('' reserved (emap.par) ndi,nri:'',2i7)')ndi,nri
c      write(*,'('' fraction'',2f5.2)')float(ima)/ndi,float(jma)/nri
c      x=float(mninp)/(jma*nd)
c      write(*,'('' mninp,  jma*nd, frac: '',2i9,f7.3)')mninp,jma*nd,x
c      x=float(mninpt)/(nr*ima)
c      write(*,'('' mninpt, nr*ima, frac: '',2i9,f7.3)')mninpt,nr*ima,x
      return
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
c force error return in debugger
      i=0
      x(1)=0.
      x(2)=0.
      xx=x(i)
      stop
c      return
      end

      subroutine errb(s,x1,x2)
      implicit real*8 (a-h,o-z)
      character*(*) s
      dimension x(2)
      write(*,'(a,1p2e12.4)')s,x1,x2
c force error return in debugger
      i=0
      x(1)=0.
      x(2)=0.
      xx=x(i)
      stop
c      return
      end

      subroutine weighd(f,f0,pha,nph,nvp,eb,iw)
c copy input weights
c f    input data
c f0   data resulting from a flat, circular, doppler map, used to calibrate weights
c      to uniform value in doppler map
c pha  phases
c nph  length f-array
c nvp  number of v-intervals
c eb   input error bars (in case they were supplied with input data)
c iw=0 internally generated weights only
c iw=1 multiply internally generated weights by weights from light curve file
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension f(*),f0(*),pha(*),eb(*)
      common /weights/ w(nd)
      common /blank/bl(npm)
      logical bl
      nil=nph*nvp
      if (nil.gt.nd) call erri('nil>nd in weighd',nil,nd)
      do i=1,nph
        do j=1,nvp
          inp=i+(j-1)*nph
          w(inp)=1
          if (f0(inp).eq.0.) w(inp)=0
          if (bl(i)) w(inp)=0
        enddo
      enddo
      if (iw.eq.1) then
c multiply by weights from error bars
        do i=1,nil
          w(i)=w(i)*f(i)/eb(i)
        enddo
      endif
c normalize to average=1
      x=0
      do i=1,nil
        x=x+w(i)
      enddo
      do i=1,nil
        w(i)=w(i)/x*nph
      enddo
      return
      end

      subroutine rldop(pht,p,pt,al,ac,flim,nim,ps,nil,nol,nit,hs,
     *           ns,ih,ipri,s,ierr)
c interface to rlfldn for doppler mapping
c RL iteration with floating default, converged to accuracy ac
c input:
c pht: input data, nil: length of input data, nol: length of solution
c p(i,j): probability p(x|xi), al: entropy weight, ac: accuracy of convergence,
c pt: transpose of p,
c flim: max acceleration parameter, nim: max no of iterations
c ns: smearing width (passed on to chims)
c ih: type of likelihood function (subroutines like, lamun)
c ipri: printout control
c out:
c ps: solution, nit: no of iterations, hs: optimum function Q.
c nit: number of iterations
c s: entropy
c ierr: 0 (OK), 1 (not converged) or 2 (not monotonic)
      implicit real*8 (a-h,o-z)
      dimension pht(*),p(*),pt(*),ps(*)
      dum=0
      idum=1
      epc=1
      call rlfldn(pht,p,pt,al,ac,flim,nim,ps,nil,nol,nit,hs,
     *           idum,idum,ns,epc,idum,dum,idum,ih,ipri,dum,s,ierr)
      return
      end

      subroutine chims(ps,chi,id1,id2,id3,ns,du1,id4,du2,id5,du3)
c default map routine to interface dopmap with emap routines
c id1-5 and du1-3 are dummies
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

      subroutine chimsi(ps,chi,id,k,nc,inc,ns,epc,ico,rdf,msu,wbkg)
c default map routine to interface dopmap with emap routines.
c since there are no sub-surfaces in dopmap, chimsi does the same thing
c as chims
      implicit real*8 (a-h,o-z)
      call chims(ps,chi,id,nc,inc,ns,epc,ico,rdf,msu,wbkg)
      return
      end

      subroutine wpar(npm,nvpm,nvm,nri,ndi)
c write a parameter file with updated values for npm,nvpm,nvm,nri,ndi
      implicit real*8 (a-h,o-z)
      open(11,file='emap.par.new')
      write(11,'(''      parameter (npm='',i3,'',nvpm='',i3,
     *           '',nvm='',i3,'')'')')npm,nvpm,nvm
      write(11,'(''      parameter (nd=npm*nvpm,nr=nvm*nvm)'')')
      write(11,'(''      parameter (nri='',i4,'',ndi='',i4,'')'')')
     *      nri,ndi
      write(11,'(''c parameters for emap routines'')')
      write(11,'(''      parameter (nf=nd,nfp=ndi,nrp=nri)'')')
      write(11,'(''      parameter (ni=nvm,nj=nvm)'')')
      write(11,'(''      parameter (nch=1,nsu=1)'')')
      close(11)
      return
      end

      function alfa(al0,ial,alf,rr,clim)
c next alfa value. last value is chosen to land rr close to clim
      implicit real*8 (a-h,o-z)
      common /alkeep/ rr2,rr1
      rr1=rr2
      rr2=rr
      if (rr.gt.sqrt(alf)*clim.or.ial.lt.3) then
        al=al0/alf**(ial-1)
      else
        dl2=log(rr1/rr2)
        dl=log(rr/clim)
        x=dl/dl2
        if (x.gt.1.) then
          al=al0/alf**(ial-1)
        else
          alfr=exp(x*log(alf))
          al=al0/alf**(ial-2)/alfr
        endif
      endif
      alfa=al
      return
      end

      subroutine chimg(x2,y2,nam,ic,nbm,jc,im,ns)
c gaussian smearing with half-width ns, by consecutive 1-d smearing in x and y
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      parameter (ms=100)
      dimension x2(ni,nj), y2(ni,nj)
      common /first/ smp(-ms:ms),ns0,ns1
      dimension t1(-ms:ni+ms),t2(-ms:nj+ms)
      if (ns0.ne.ns) then
c (re)make smearing matrix smp if value of ns is new.
        w=1.5
c actual smearing width used is w times 1/e width
        ns1=w*ns
        if (ns1.gt.ms) call errb('ns1>ms in chimg',ns1,ms)
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
c RL iteration with floating default, converged to accuracy ac
c input:
c pht: input data, nil: length of input data, nol: length of solution
c p(i,j): probability p(x|xi), al: entropy weight, ac: accuracy of convergence,
c flim: max acceleration parameter, nim: max no of iterations
c nil: no of phases, nol: size of image, nc: no of components of default,
c epc: the weights of each of these components, inc: indices for the types of
c component.
c out:
c ps: solution, nit: no of iterations, hs: optimum function Q.
c ih: type of likelihood function (subroutines like, lamun)
c ierr: 0 (OK), 1 (not converged) or 2 (not monotonic)
c ipri: printout control
      implicit real*8 (a-h,o-z)
      real*4 dpx
      include 'emap.par'
      dimension pht(*),ps(*),inc(nch,nsu),ns(nch,nsu),epc(*),ico(*)
      real*4 p(nfp,nr),pt(nrp,nf)
      dimension ph(nf),dha(nr),dsa(nr),chi(nr,nch),tem(nr),
     *          ex(nr),dps(nr),dps0(nr)
      common /difpix/ dpx(nr,4)
      ierr=0
      del=1
      nit=0
      hs=0
c sum of light curve
      spt=0
      do i=1,nil
        spt=spt+pht(i)
      enddo
      do while (del.gt.ac.and.nit.lt.nim)
        nit=nit+1
c phi (forward projection step)
        call fproj(pt,nil,nol,ps,ph)
c magic statement: On HP's, code crashes with bus error if this line is removed:
        if (nit.lt.0) write(*,*)1
        if (al.le.0.) then
c measure convergence by diff between ph and pht
          av=0
          sm=0
          do i=1,nil
            av=av+pht(i)
            sm=sm+(pht(i)-ph(i))**2
          enddo
          rm=sqrt(sm/nil)
          av=av/nil
          del=rm/av
        endif
c H and dH/dpsi
        call like(pht,ph,p,nil,nol,ih,h,dha)
c sum of psi*dH/dpsi and sum of psi
        su=0
        sup=0
        do j=1,nol
          su=su+ps(j)*dha(j)
          sup=sup+ps(j)
        enddo
c delta-H (in same array as dha)
        do j=1,nol
          dha(j)=dha(j)+(-su/sup+spt/sup-1)
        enddo
c entropy
        if (al.gt.0.) call entr(ps,nil,nol,s,dsa,chi,tem,ex,nc,inc,ns,
     *       epc,ico,rdf,msu,wbkg)
c sum of psi*dS/dpsi
        su=0
        do k=1,nol
          su=su+ps(k)*dsa(k)
        enddo
c Delta-S (in same array as dsa)
        do j=1,nol
          dsa(j)=dsa(j)-su/sup
        enddo
c save old value of Q
        hs0=hs
c current value of Q
        hs=h+al*s
c warning if not monotonic
        if (nit.ne.1.and.hs.lt.hs0) ierr=2
c new dpsi
        do j=1,nol
          dps(j)=ps(j)*(dha(j)+al*dsa(j))
        enddo
        if (nit.eq.1) then
          do j=1,nol
            dps0(j)=dps(j)
          enddo
        endif

c determine optimum direction and factor for update, from dps0 and dps.
c (generalization of Lucy's acceleration)
c derivatives of Q w/r to lambda and mu at current psi
        call lamun(pht,ps,dps0,dps,pt,ph,chi,tem,ex,al,nil,nol,
     *             fla,fmu,nc,inc,ns,epc,ico,rdf,msu,ih,wbkg)
c length of multiplier
        amu=sqrt(fla**2+fmu**2)
c direction (normalized) of optimum
        fla=fla/amu
        fmu=fmu/amu
c find max of amu
        flm=0
        jlm=0
        do j=1,nol
c value of prospective update
          xx=fla*dps0(j)+fmu*dps(j)
          if (xx.lt.0) then
c positivity not maintained at this value of amu, get max
            fl=-xx/ps(j)
            if (fl.gt.flm) then
              flm=fl
              jlm=j
            endif
          endif
        enddo
        if (flm.gt.0.) then
c positivity restriction limits amu
c safety factor 0.9
          flm=dmin1(flim,0.9/flm)
        else
c take overall limit
          flm=flim
        endif
        amu=dmin1(amu,flm)
c check if amu has been small for awhile
        call zcheck(amu,nit,20,0.01d0,ierr)
        if (ierr.eq.4) return
c having found amu, update psi
        do j=1,nol
          up=fmu*dps(j)+fla*dps0(j)
          ps(j)=ps(j)+amu*up
c save direction of old dps
          dps0(j)=up
        enddo

        if (al.gt.0.) then
          del=0
          do j=1,nol
c measure of convergence (Lucy's eq 17)
            dj=abs(dha(j)+al*dsa(j))/(abs(dha(j))+al*abs(dsa(j)))
            del=del+dj**2
          enddo
          del=sqrt(del/nol)
        endif
c store change
        do j=1,nol
          dpx(j,mod(nit,4)+1)=dha(j)+al*dsa(j)
        enddo
c printout
        if (nit.eq.1)  then
          write(10,'(''  it        H+alfa*S        delta'')')
          if (ipri.gt.1)
     *      write(*,'('' it      H+alfa*S        delta'')')
        endif
        if (nit.ne.1.and.hs.lt.hs0) then
c line summary to activity log file
          write(10,'(i4,1pe20.12,e11.2,'' **'')')nit,hs,del
          if (ipri.gt.1)
     *       write(*,'(i4,1pe20.12,e11.2,'' **'')')nit,hs,del
        else
          write(10,'(i4,1pe20.12,e11.2)')nit,hs,del
          if (ipri.gt.1)
     *       write(*,'(i4,1pe20.12,e11.2)')nit,hs,del
        endif
c  renormalize (mostly to strip accumulating numerical noise)
        su=0
        do j=1,nol
          su=su+ps(j)
        enddo
        do j=1,nol
          ps(j)=ps(j)*spt/su
        enddo
      enddo
      if (nit.ge.nim) then
c not converged at this value of ac
        hs=0
c return last value of del
        ac=del
c error flag
        ierr=ierr+1
      endif
      return
      end

      subroutine like(pht,ph,p,nil,nol,ih,h,dha)
c likelihood function h and its derivative w/r psi, dha
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension pht(*),ph(*),dha(*)
      real*4 p(nfp,nr)
      dimension phh(nf)
      common /weights/ w(nf)
      h=0
      if (ih.eq.0) then
c log likelihood
        do i=1,nil
          h=h+w(i)*(pht(i)*log(ph(i))-ph(i))
        enddo
      else
c rms
        do i=1,nil
          h=h-w(i)*(pht(i)-ph(i))**2
        enddo
c normalizing factor
        h=h*nil
      endif
c Delta-H
      if (ih.eq.0) then
c log likelihood
        do i=1,nil
          phh(i)=w(i)*(pht(i)/ph(i)-1)
        enddo
        call bproj(p,nil,nol,phh,dha)
      else
c rms
        do i=1,nil
          phh(i)=2*nil*w(i)*(pht(i)-ph(i))
        enddo
        call bproj(p,nil,nol,phh,dha)
      endif
      return
      end

      subroutine zcheck(amu,nit,nz,amc,ierr)
c set ierr=4 if past nz amu values were less than amc
      implicit real*8 (a-h,o-z)
      common /zsave/ am(30)
      if (nz.gt.30) call errb('nz>30 in zcheck',nz,30)
      am(mod(nit,nz)+1)=amu
      if (nit.ge.nz) then
        ama=0
        do i=1,nz
          if(am(i).gt.ama) ama=am(i)
        enddo
        if (ama.lt.amc) ierr=4
      endif
      return
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
      if (nc.gt.nch) call errb('nc>nch in entr',nc,nch)
c make chi's
      call chims(ps,chi,1,nc,inc,ns,epc,ico,rdf,msu,wbkg)
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
          call chimsi(x,y,-1,k,nc,inc,ns,epc,ico,rdf,msu,wbkg)
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
      implicit real*8 (a-h,o-z)
      include 'emap.par'
      dimension pht(*),ps(*),dps1(*),dps2(*),ph(*),chi(nr,nch),
     *          tem(*),ex(*),inc(nch,nsu),ns(nch,nsu),epc(*),ico(*)
      real*4 pt(nrp,nf)
      dimension dph1(nf),dph2(nf),dch1(nr,nch),dch2(nr,nch)
      common /weights/ w(nf)
c changes in phi and chi corresponding to dps0, dps
      call fproj(pt,nil,nol,dps1,dph1)
      call fproj(pt,nil,nol,dps2,dph2)
      if (al.gt.0.) then
        call chims(dps1,dch1,1,nc,inc,ns,epc,ico,rdf,msu,wbkg)
        call chims(dps2,dch2,1,nc,inc,ns,epc,ico,rdf,msu,wbkg)
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
c likelihood term
      if (ih.eq.0) then
c log likelihood
        do i=1,nil
          phw=pht(i)*w(i)
          dpp=pht(i)/ph(i)-1
          xl=dph1(i)/ph(i)
          xm=dph2(i)/ph(i)
          hl=hl+w(i)*dph1(i)*dpp
          hll=hll-xl**2*phw
          hlm=hlm-xl*xm*phw
          hm=hm+w(i)*dph2(i)*dpp
          hmm=hmm-xm**2*phw
        enddo
      else
c rms likelihood
        do i=1,nil
          x=pht(i)-ph(i)
          xl=dph1(i)
          xm=dph2(i)
          hl=hl+2*w(i)*x*xl
          hll=hll-2*w(i)*xl**2
          hlm=hlm-2*w(i)*xl*xm
          hm=hm+2*w(i)*x*xm
          hmm=hmm-2*w(i)*xm**2
        enddo
        hl=hl*nil
        hm=hm*nil
        hll=hll*nil
        hmm=hmm*nil
        hlm=hlm*nil
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
c calculate lambda and mu
      dd=qmm*qll-qlm**2
      dr=abs(dd)/(qll**2+qmm**2)
c if det of 2X2 system too small, use only mu
      if (dr.lt.1e-3) then
        fla=0
        fmu=-qm/qmm
      else
        fmu=(qlm*ql-qll*qm)/dd
        fla=(qlm*qm-qmm*ql)/dd
      endif
      return
      end

      function rrms(f1,eb,iw,f0,n)
c rms of difference between f1 and f0 (iw=0) or
c rms distance between f1 and f0 normalized by the error bars eb (iw=1)
      implicit real*8 (a-h,o-z)
      dimension f1(*),eb(*),f0(*)
      s=0
      if (iw.eq.0) then
        do i=1,n
          s=s+(f1(i)-f0(i))**2
        enddo
      else
        do i=1,n
          s=s+((f1(i)-f0(i))/eb(i))**2
        enddo
      endif
      rrms=sqrt(s/(n-1))
      return
      end

      function rmsn(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common/folds/nph,nvp
      s=0
      do i=2,nph-1
        do j=2,nvp-1
          i1=i+(j-2)*nph
          i2=i1+nph
          i3=i2+nph
          s0=x(i1-1)+x(i1)+x(i1+1)+x(i2-1)+x(i2+1)+x(i3-1)+x(i3)+x(i3+1)
          s1=y(i1-1)+y(i1)+y(i1+1)+y(i2-1)+y(i2+1)+y(i3-1)+y(i3)+y(i3+1)
          s=s+((s1-s0)/8-(y(i2)-x(i2)))**2
        enddo
      enddo
      s=s/(nph-2)/(nvp-2)
      rmsn=sqrt(s)
      return
      end
