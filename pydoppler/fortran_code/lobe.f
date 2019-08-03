      program stream
c aux program for idl plots of dopmaps
c calculates Roche lobes and integrates path of stream from L1
      parameter(nmax=10000)
      complex z,w,z1,z2,dz,dw,wk,no,wr,wc,w0,wm1,wi,
     *       wout(nmax),wkout(nmax)
      dimension xout(nmax),yout(nmax),rout(nmax)
      open(3,file='lobe.out')
      open(4,file='lobe.in')
      read(4,*)qm
      if (abs(qm-1).lt.1e-4) qm=1e-4
      close(4)
      rd=0.1
      if (qm.le.0.) then 
        write(3,'(''qm<0: '',f8.5)')qm
        close(3)
        stop
      endif
      rl1=rlq1(qm)
c      write(*,'('' rlq1'',f8.5)')rl1
      write(3,'(f8.5)')qm
      call lobes(qm,rl1)
c center of mass rel to M1
      cm=qm/(1+qm)
c coordinates of M1 and M2
      z1=-cm
      z2=1-cm
      wm1=conjg(cmplx(0.,-cm))
c start at L1-eps with v=0
      eps=1e-3
      z=cmplx(rl1-cm-eps,0.)
      w=0
c check that we really are near L1
      call eqmot(z,w,zp,wp,z1,z2,qm)
      write(*,'('' t=0: z,w'',1p4e11.3)')zp,wp
      t=0
      dt=1e-4 
      isa=0
      it=0
      r=1
      ist=0
      ph=0.
      phmax=6
c      do while(it.lt.nmax.and.ph.lt.phmax.and.ist.eq.0)
      do while(it.lt.nmax.and.ph.lt.phmax)
        it=it+1
        call intrk(z,w,dt,dz,dw,z1,z2,qm) 
        z=z+dz
        w=w+dw 
        t=t+dt
        if (abs(dz)/abs(z).gt.0.02) dt=dt/2
        if (abs(dz)/abs(z).lt.0.005) dt=2*dt
        dph=-aimag(z*conjg(z-dz))/abs(z)/abs(z-dz)
        ph=ph+dph
c velocity in inertial frame
c change by Guillaume
        wi=w+cmplx(0,1.)*z
c unit vector normal to kepler orbit
        rold=r
        r=abs(z-z1)
        if (ist.eq.0.and.rold.lt.r) then
          ist=1
          rmin=rold
        endif
c kepler velocity of circular orbit in potential of M1, rel. to M1
        vk=1/sqrt(r*(1+qm))
c unit vector in r
        no=conjg(z-z1)/r
        wk=-vk*no*cmplx(0.,1.)
c same but rel. to cm, this is velocity in inertial frame
        wk=wk+wm1
c velocity normal to disk edge, in rotating frame
        dot=no*w
c velocity parallel to disk edge
        par=aimag(no*w)
c reflected velocity
        wr=w-2*dot*no
c        write(*,'(f8.4,1p9e11.3)')t,z,w,wk,wr,r
        xout(it)=z+cm
        yout(it)=-aimag(z)
        rout(it)=sqrt(xout(it)**2+yout(it)**2)
c change by Guillaume
        wout(it)=wi								
        wkout(it)=conjg(wk)
        if (it.gt.1) then 
          x=xout(it)
          y=yout(it)
          phi=atan(y/x)
          if (rout(it).lt.rd.and.rout(it-1).gt.rd) then
c            write(*,'('' r,x,y,phi,vs,vk,dot,par'',8f8.3)')
c     *         rout(it),x,y,phi,real(w),vk,dot,par
c            write(*,'('' w,no'',4f8.3)')w,no
            x=xout(it-1)
            y=yout(it-1)
            phi=atan(y/x)
c            write(*,'('' r,x,y,phi'',4f8.3)')rout(it-1),x,y,phi
          endif
        endif
        if (isa.eq.0.and.yout(it).lt.0) then
          isa=1
          ra=abs(z-z1)
          wc=conjg(w)+cmplx(0.,1.)*conjg(z-z1)
          ang=abs(aimag((z-z1)*conjg(wc)))
        endif
      enddo
      write(3,'(i5)')it
      write(3,'(1p5e15.7)')(xout(i),i=1,it)
      write(3,'(1p5e15.7)')(yout(i),i=1,it)
      write(3,'(1p6e12.4)')(real(wout(i)),i=1,it)
      write(3,'(1p6e12.4)')(aimag(wout(i)),i=1,it)
      write(3,'(1p6e12.4)')(real(wkout(i)),i=1,it)
      write(3,'(1p6e12.4)')(aimag(wkout(i)),i=1,it)
      rc=ang**2*(1+qm)
      z=rl1-cm
      w0=cmplx(0.,1.)*conjg(z-z1)
      ang0=abs(aimag((z-z1)*conjg(w0)))
      rc0=ang0**2*(1+qm)
      write(3,'(3f8.5,'' ra, rc '')')ra,rc,rc0
c      write(*,'(3f8.5,'' ra, rc '')')ra,rc,rc0
c      write(*,'(2f8.4,'' q, rmin'')')qm,rmin
      end

      subroutine intrk(z,w,dt,dz,dw,z1,z2,qm)
      complex z,w,dz,dw,z1,z2
      complex zp,wp,hz0,hz1,hz2,hz3,hw0,hw1,hw2,hw3,zx,wx
      zx=z
      wx=w
      call eqmot(zx,wx,zp,wp,z1,z2,qm)
      hz0=zp*dt
      hw0=wp*dt
      zx=z+hz0/2
      wx=w+hw0/2
      call eqmot(zx,wx,zp,wp,z1,z2,qm)
      hz1=zp*dt
      hw1=wp*dt
      zx=z+hz1/2
      wx=w+hw1/2
      call eqmot(zx,wx,zp,wp,z1,z2,qm)
      hz2=zp*dt
      hw2=wp*dt
      zx=z+hz2
      wx=w+hw2
      call eqmot(zx,wx,zp,wp,z1,z2,qm)
      hz3=zp*dt
      hw3=wp*dt
      dz=(hz0+2*hz1+2*hz2+hz3)/6
      dw=(hw0+2*hw1+2*hw2+hw3)/6
      return
      end      

      subroutine eqmot(z,w,zp,wp,z1,z2,qm)
      complex z,w,zp,wp,z1,z2
      complex zr1,zr2
      zr1=z-z1
      zr2=z-z2
c change by Guillaume : - sign in Coriolis
      wp=-(qm*zr2/(abs(zr2))**3+zr1/(abs(zr1))**3)/(1+qm)-
     *    cmplx(0.,2.)*w+z
      zp=w
      return
      end
    
      function rlq1(q)
c radius (from primary) of L1
      if (abs(1-q).lt.1e-4) then
        rlq=0.5
        return
      endif
      rl=0
      rn=1-q
 2    if (abs(rl/rn-1).gt.1e-4)then
        rl=rn
        f=q/(1-rl)**2-1/rl**2+(1+q)*rl-q
        fa=2*q/(1-rl)**3+2/rl**3+(1+q)
        rn=rl-f/fa
        goto 2
      endif
      rlq1=rn
      return
      end

      function pot(q,x,y,z,pr)
c Roche potential. coordinates centered on M2,
c z along rotation axis, x toward M1
c pr is gradient in radius from M2
c first transform to polar coordinates w/r rotation axis
      r=sqrt(x*x+y*y+z*z)
      if (r.eq.0) stop 'r=0 in pot'
      rh=sqrt(x*x+y*y)
      st=rh/r
      if (rh.eq.0) then
        cf=1
      else
        cf=x/rh
      endif
      r2=1/(1+q)
      r1=sqrt(1+r**2-2*r*cf*st)
      pot=-1/r-1/q/r1-0.5*(1/q+1)*(r2**2+(r*st)**2-2*r2*r*cf*st)
      pr=1/r**2+1/q/(r1**3)*(r-cf*st)-0.5*(1/q+1)*2*(r*st*st-r2*cf*st)
      return
      end

      subroutine lobes(q,rs)
      include 'emap.par'
      dimension r(ni,nj),ch(ni),ps(nj),x(ni),y(ni)
      nc=ni
      np=nj
      call surf(q,rs,nc,np,r,ch,ps)
      write(3,'(2i5,'' ni,nj'')')ni,nj
      j=1
      do i=1,nc
        x(i)=1-r(i,j)*cos(ch(i))
        y(i)=-r(i,j)*sin(ch(i))
      enddo
      write(3,'(8f9.5)')(x(i),i=1,nc),(x(i),i=nc,1,-1)
      write(3,'(8f9.5)')(y(i),i=1,nc),(-y(i),i=nc,1,-1)
      call surf(1/q,1-rs,nc,np,r,ch,ps)
      j=1
      do i=1,nc
        x(i)=r(i,j)*cos(ch(i))
        y(i)=r(i,j)*sin(ch(i))
      enddo
      write(3,'(8f9.5)')(x(i),i=1,nc),(x(i),i=nc,1,-1)
      write(3,'(8f9.5)')(y(i),i=1,nc),(-y(i),i=nc,1,-1)
      return
      end

      subroutine surf(q,rs,nc,np,r,ch,ps)
c Roche surface around M2, coordinates on surface are ch, ps.
c ch: polar angle from direction to M1; ps: corresponding azimuth, counting
c from orbital plane.
c q:mass ratio, rs: radius of surface at point facing M1
c nc, np: number of chi's, psi's.
c output:
c r(nf,nt): radius. ch, ps: chi and psi arrays
      include 'emap.par'
      dimension r(ni,nj),ch(*),ps(*)
      if (nc.gt.ni) stop 'nc>ni'
      if (np.gt.nj) stop 'np>nj'
      pi=2*asin(1.)
      dc=pi/nc
      ch(1)=0
      do i=1,nc 
        ch(i)=float((i-1))*pi/(nc-1)
      enddo
      ps(1)=0
      do j=1,np
        ps(j)=float((j-1))*2*pi/np
      enddo
      rs1=1-rs
      fs=pot(q,rs1,0.,0.,pr)
c max no of iterations
      im=20
      do i=1,np
        cp=cos(ps(i))
        sp=sin(ps(i))
        rx=(1-dc)*rs1
        r(1,i)=rs1
        do k=2,nc
          x=cos(ch(k))
          sc=sin(ch(k))
          y=sc*cp
          z=sc*sp
          j=0
          f=1
          do while (j.lt.im.and.abs(f-fs).gt.1e-4.or.j.eq.0.)
            j=j+1
            r1=rx
            f=pot(q,r1*x,r1*y,r1*z,pr)
            rx=r1-(f-fs)/pr
            if (rx.gt.rs1) rx=rs1
          enddo  
          if (j.ge.im) then
            write(*,'('' no conv in surf; k,i,ch,ps'',2i4,2f7.3)')
     *           k,i,ch(k),ps(i)
            stop
          endif
          r(k,i)=rx  
        enddo 
      enddo   
      return
      end
