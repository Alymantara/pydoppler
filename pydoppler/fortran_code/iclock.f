      subroutine clock(t)
      real*8 t
c ibm-rs clock function
      t=mclock()/100.
      return
      end
