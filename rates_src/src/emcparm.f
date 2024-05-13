      SUBROUTINE EMCPARM(taremc,x,emc)

      implicit none
      integer*4 taremc
      real*8 x,emc

CCCCCCCC    H2     CCCCCCCCC

      if(taremc.EQ.1) emc = 1.

CCCCCCCC    D2     CCCCCCCCC

      if(taremc.EQ.2) emc = 1.

CCCCCCCC    D2     CCCCCCCCC

      if(taremc.EQ.3) emc = 1.

CCCCCCCC    He     CCCCCCCCC

      if(taremc.EQ.4) emc = 1.

CCCCCCCC    C     CCCCCCCCC

      if(taremc.EQ.5) emc = (0.926-0.400*x-0.0987*exp(-27.714*x)+ 
     &                 0.257*x**0.247)  

CCCCCCCC    Al     CCCCCCCCC

      if(taremc.EQ.6) emc = (0.825-0.46*x-0.19*exp(-21.8*x)+
     &                    (0.34*x**(-4.91))*x**5.0)

CCCCCCCC    Cu     CCCCCCCCC

      if(taremc.EQ.7) emc = (1.026-0.56*x-0.34*exp(-45.7*x)+(0.26*
     &                   x**(-4.41))*x**5.0)
      

CCCCCCCC    Fe     CCCCCCCCC

      if(taremc.EQ.8) emc = 1.


CCCCCCCC    Au     CCCCCCCCC

      if(taremc.EQ.9) emc = (0.970-1.433*x-0.334*exp(-54.53*x)+
     &                  1.074*x**0.711)


      return
      end
  
