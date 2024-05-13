         SUBROUTINE i_d2_model(Q2,W2,SIG_NRES,SIG_RES1,SIG_RES2,
     >                    sigroper,goroper)
********************************************************************************
*
* This subroutine calculates model cross sections for deuterium in resonance region. 
* Cross section is returned in nanobarns. 
* For all the resonance regions we use nonrelativistc
* Breit-Wigner shapes (L.Stuart). 
* This subroutine is based on h2model.f from SLAC
* 
*
* E        = Incident electron energy in GeV.
* EP       = Scattered electron energy in GeV.
* TH       = Scattered electron angle in degrees.
********************************************************************************
      IMPLICIT NONE

      INCLUDE 'imodel.cmn'
*
*      REAL*8 XVAL(37)
* 
      logical goroper
      INTEGER I
*
* modified by I.N. 04/09/97
*

*      LOGICAL FIRST
      REAL*8 W2, Q2
      REAL*8 SIG_NRES(2),SIG_RES1(2),SIG_RES2(2),sigroper(2), 
     >     WR, KLR, KRCM, EPIRCM, PPIRCM, GG, GPI, DEN, 
     >     SIGDEL, W, KCM, K, EPICM, PPICM, WDIF, GAM, 
     >     GR, RERRMAT(25,25), RSIGERR(25),ERRMAT(22,22), SIGERR(22), 
     >     SIGTOTER, ERRCHECK

      REAL*8 PI, ALPHA, AM,
     >     MPPI, MDELTA, MPI, 
     >     GAM1WID, GAM2WID, MASS1,
     >     ROPERWID, MASSROPER, 
     >     MASS2, DELWID, FIT(30), SQRTWDIF,
     >     XR, MQDEP

      integer jj
      real*8 kr,meta
      real*8 mrho,mn,mn2,am2
      real*8 gamma_gamma,gamma_pi,denom_pi
      real*8 sigma_delta,sigma_second,sigma_third
*
      real*8 mass(3),width(3),qdep(3)
      real*8 coeff(3,0:3),nres_coeff(0:4,3)
* Define vectors : MASS(3), WIDTH(3), which will contain the masses and widths of
* the three resonances: 1=DELTA, 2=second res reg, 3=third res reg
*
      DATA PI/3.14159265/
      DATA ALPHA/0.00729735/
      DATA AM/.938259/
      DATA MPI/0.13957/
      DATA MPPI/1.07783/ 
      data Meta/0.54745/ 
      data Mrho/0.7681/
*      DATA XR/0.18/
	data Mn/0.938/
      Am2 = Am*Am
*
* The fit parameters are called  XVAL 
* Assign them to masses, widths, etc.
*
      xr=XVAL(1)
      do i=1,3
         mass(i)  = XVAL(i+1)
*         write(*,*)xval(i)
      enddo
      do i=1,3
         width(i)  = XVAL(i+4)
*         write(*,*)xval(i+3)
      enddo
      do i=1,3
         qdep(i)   = XVAL(i+7)
*         write(*,*)xval(i+6)
      enddo     
      k=11
      do i=1,3
         do j=0,3
            coeff(i,j) = XVAL(k)
*            write(*,*)xval(k)
            k=k+1
         enddo
      enddo
      do i=0,4
         do j=1,3
           nres_coeff(i,j) = XVAL(k)
*          write(*,*)xval(k)
           k=k+1
         enddo
      enddo
*      write(*,*)'Done assigning parameters'
*
      W = SQRT(W2) 
      WDIF = MAX(0.0001,W - MPPI)
****************************************************************
* This part is not used for deuterium.
****************************************************************
* K equivalent real photon energy needed to excite resonance W
* K is in Lab frame, KCM is in CM frame
* 
         K = (W2 - Am2)/2./Am
         KCM = (W2 - Am2)/2./W
         KR = (Mass(1)*Mass(1)-Am2)/2./Am
         KRCM = (Mass(1)*Mass(1) - Am2)/2./Mass(1)
* Decay mode for Delta: N,pi:
********
      EPICM = 0.5*( W2 + Mpi*Mpi - Am2 )/W

      PPICM = SQRT(MAX(0.0,(EPICM*EPICM - MPI*MPI)))
      EPIRCM = 0.5*(Mass(1)*Mass(1) + Mpi*Mpi - Am2 )/Mass(1)   
      PPIRCM = SQRT(MAX(0.0,(EPIRCM*EPIRCM - MPI*MPI)))
      
********
* Now calculate the partial widths:
*  gamma_pi(i)
*      write(*,*)'mass(1),w2 ',mass(1),' ' ,w2
*      write(*,*)'width(1) ',width(1)
*      write(*,*)'PPICM, PPIRCM ',PPICM,' ',PPIRCM
*      write(*,*)'XR ' ,xr
      gamma_pi = width(1)*(PPICM/PPIRCM)**3*(PPIRCM*PPIRCM+XR*XR)/
     >          (PPICM*PPICM+XR*XR)
*
*      write(*,*)'iuhu 2'
      gamma_gamma = width(1)*(KCM/KRCM)**2*(KRCM*KRCM+XR*XR)/
     >          (KCM*KCM+XR*XR)

      denom_pi = (W2 - Mass(1)*Mass(1))**2 + (Mass(1)*gamma_pi)**2
*      write(*,*)'iuhu 3'
***********************************************************************
***********************************************************************
* Now calculate cross sections for Delta:
      sigma_delta = width(1)/((W-Mass(1)*(1. + Q2*qdep(1)))**2 
     >               + 0.25*width(1)*width(1))
* For the second and third resonance regions:
      sigma_second = width(2)/((W-Mass(2)*(1. + Q2*qdep(2)))**2 
     >               + 0.25*width(2)*width(2))

*      write(*,*)'width(2)',' ',width(2)
*      write(*,*)'mass(2)',' ',mass(2)

*      write(*,*)'sigma_second',' ',sigma_second

      sigma_third  = width(3)/((W-Mass(3)*(1. + Q2*qdep(3)))**2
     >               + 0.25*width(3)*width(3))

*      write(*,*)'sigma_third',' ',sigma_third

*      write(*,*)'width(3)',' ',width(3)
*      write(*,*)'mass(3)',' ',mass(3)

*      write(*,*)'iuhu 5'
*
* Put in the Q2 dependence of AH(Q2):
* AH(Q2) = coeff(i,j)*Q2**j
*
      SIG_RES1(1) = 0.0
      SIG_RES2(1) = 0.0
      SIG_NRES(1) = 0.0
      do j=0,3  
         sig_res1(1) = sig_res1(1)+ sigma_delta*coeff(1,j)*Q2**j
         sig_res2(1) = sig_res2(1)+
     >               sigma_second*coeff(2,j)*Q2**j+
     >               sigma_third *coeff(3,j)*Q2**j
      enddo   
*      write(*,*)'sigma_d,sigma_23',' ',sig_res1(1),' ',sig_res2(1)
      do j=0,4
         do jj=1,3
            sig_nres(1) = sig_nres(1)+
     >          nres_coeff(j,jj)*q2**j*sqrt(wdif**(2*jj-1))
         enddo
      enddo     
*
      RETURN
      END




