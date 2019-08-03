      subroutine clock(t)
      real*8 t
      real*4 ta(2)
      x=dtime(ta)
      t=t+ta(1)
      return
      end
