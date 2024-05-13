      program WISER
! get p, pi, or K cross sections from Wiser fit
! valid for 5<E<19 GeV and 5<theta<50 degrees
! from Steve Rock

      implicit none
      real E0,P,THDEG,RL,SIGPIP,SIGPIM,SIGKP,SIGKM,SIGPP,SIGPB
      real LHmult,LDmult,dSigma,dP
      real dsig,I,Ne,Na,Np
      real rpip,rpim,rkp,rkm,rpp,rpb,rptot,rmtot
      real xsigpip(10),xsigpim(10),xsigkp(10),xsigkm(10),xsigpp(10)
      real xsigpb(10),xp(10),xj
      integer ispec,j

 5    write(6,*)' 1=HMS; 2=SHMS '
      read(5,*,err=900) ispec

      if (ispec.eq.1) then      ! HMS
         dSigma=6
         dP = 0.20              ! +/-10%
         write(6,*)' rates will be for HMS'
      elseif (ispec.eq.2) then  ! SHMS
         dSigma=4
         dP = 0.40              ! +/-20%
         write(6,*)' rates will be for SHMS'
      else
         write(6,*)' invalid spectrometer choice'
      endif

 10   write(6,*)' Input: '
      write(6,*)'   Electron beam energy (5-19 GeV): '
      read(5,*,err=900)E0
      write(6,*)'   Scattered particle momentum (GeV/c): '
      read(5,*)P
      write(6,*)'   Scattered angle (5-50 deg): '
      read(5,*)THDEG
      write(6,*)'   Radiation length of target (%): '
      read(5,*)RL

      write(6,100)E0,P,THDEG,RL
 100  format(/,' E0=',f6.3,'  P=',f6.3,'  TH=',f6.3,'  RL=',f5.3)
      write(6,101)
 101  format(4x,'pi+',9x,'pi-',10x,'K+',10x,'K-',11x,'p',8x,'pbar',
     *     3x,'nb/GeV-sr')

c calculate for central momentum only
! pi+
c      CALL WISER_ALL_SIG (E0,P,THDEG,RL,1,SIGPIP)
! pi-
c      CALL WISER_ALL_SIG (E0,P,THDEG,RL,2,SIGPIM)
! K+
c      CALL WISER_ALL_SIG (E0,P,THDEG,RL,3,SIGKP)
! K-
c      CALL WISER_ALL_SIG (E0,P,THDEG,RL,4,SIGKM)
! p
c      CALL WISER_ALL_SIG (E0,P,THDEG,RL,5,SIGPP)
! p-bar
c      CALL WISER_ALL_SIG (E0,P,THDEG,RL,6,SIGPB)

c gh - 18.10.12
c subdivide spectrometer focal plane into 10 pieces to get a more
c accurate average cross section over the focal plane

      sigpip=0.
      sigpim=0.
      sigkp=0.
      sigkm=0.
      sigpp=0.
      sigpb=0.
      
      do j=1,10
         xj=float(j)
         xp(j)=p*(1.+dp/10.*(5.5-xj))

! pi+
         CALL WISER_ALL_SIG (E0,xP(j),THDEG,RL,1,xSIGPIP(j))
! pi-
         CALL WISER_ALL_SIG (E0,xP(j),THDEG,RL,2,xSIGPIM(j))
! K+
         CALL WISER_ALL_SIG (E0,xP(j),THDEG,RL,3,xSIGKP(j))
! K-
         CALL WISER_ALL_SIG (E0,xP(j),THDEG,RL,4,xSIGKM(j))
! p
         CALL WISER_ALL_SIG (E0,xP(j),THDEG,RL,5,xSIGPP(j))
! p-bar
         CALL WISER_ALL_SIG (E0,xP(j),THDEG,RL,6,xSIGPB(j))

         sigpip=sigpip+xsigpip(j)
         sigpim=sigpim+xsigpim(j)
         sigkp=sigkp+xsigkp(j)
         sigkm=sigkm+xsigkm(j)
         sigpp=sigpp+xsigpp(j)
         sigpb=sigpb+xsigpb(j)
      enddo

      sigpip=sigpip/10.
      sigpim=sigpim/10.
      sigkp=sigkp/10.
      sigkm=sigkm/10.
      sigpp=sigpp/10.
      sigpb=sigpb/10.

c      do j=1,10
c         write(6,'(i3,7e12.4)')j,xp(j),xSIGPIP(j),xSIGPIM(j),xSIGKP(j),
c     a     xSIGKM(j),xSIGPP(j),xSIGPB(j)
c      enddo
      
      write(6,'(6e12.4)')SIGPIP,SIGPIM,SIGKP,SIGKM,SIGPP,SIGPB
      
      write(6,*)' Rates '
!============================================================================
! gh 10.12.10
!     Konrad Aniol and I confirm that LHmult=98 is 90ua on 4cm LH2 tgt,
!     if 0.2DeltaP/P momentum bite and 5msr solid angle for the SHMS are
!      assumed
!============================================================================

!============================================================================
! Samip 25.05.16
!     LHmult for 70uA on 10cm LH2 tgt is 190 similar to the calculation above
!      for 90 uA on 4cm tgt
!     Also, the most updated solid angle acceptance is 4msr and momentum
!      acceptance is -10% to 22%. From Bill's suggestion, using 15% times 2
!============================================================================

!============================================================================
! Samip 25.10.17
! LHmult for 20uA on 10cm LH2 tgt is ~55.
!============================================================================

!============================================================================
!Samip 25.10.17
! LHmult for 50uA on 10cm LH2 tgt is ~136
!============================================================================

!============================================================================
!Samip 01.11.17
! LHmult for 5 uA on 3 % radiation length Carbon target is ~24
!============================================================================

      LHmult=190
      rpip = SIGPIP*LHmult*dSigma*dP*P
      rpim = SIGPIM*LHmult*dSigma*dP*P
      rkp  = SIGKP*LHmult*dSigma*dP*P
      rkm  = SIGKM*LHmult*dSigma*dP*P
      rpp  = SIGPP*LHmult*dSigma*dP*P
      rpb  = SIGPB*LHmult*dSigma*dP*P
      rptot = rpip+rkp+rpp
      rmtot = rpim+rkm+rpb

      write(6,200)rpip,rpim,rkp,rkm,rpp,rpb
      write(6,201)rptot,rmtot
 200  format('H',8e12.4)
 201  format('T',2e12.4)

!     LDmult=115
!      write(6,201)SIGPIP*P*LDmult,SIGPIM*P*LDmult,SIGKP*P*LDmult,
!      *	SIGKM*P*LDmult,SIGPP*P*LDmult,SIGPB*P*LDmult
! 201  format('D',6e12.4)
      
C Method from Dave Gaskell is more transparent
C step 1 convert nb to cm2  *1E-09 barn/nb * 1E-28 m**2/barn *(100cm)**2/1 m**2
C step 2 convert current to number of electrons
      
      I  = 70.0              ! in uA
      Ne = I*1.0E-6/1.6E-19  ! 1/s
      Na = 6.022E23
      Np = 10.0*0.0708*Na/1.008 ! 1/cm**2
      dsig=dSigma/1000.         ! convert msr to sr

      rpip = sigpip*1.0E-9*1.0E-28*100.0*100.0*Ne*Np*dP*P*dsig
      rkp  =  sigkp*1.0E-9*1.0E-28*100.0*100.0*Ne*Np*dP*P*dsig
      rpp  =  sigpp*1.0E-9*1.0E-28*100.0*100.0*Ne*Np*dP*P*dsig
      
      rpim = sigpim*1.0E-9*1.0E-28*100.0*100.0*Ne*Np*dP*P*dsig
      rkm  =  sigkm*1.0E-9*1.0E-28*100.0*100.0*Ne*Np*dP*P*dsig
      rpb  =  sigpb*1.0E-9*1.0E-28*100.0*100.0*Ne*Np*dP*P*dsig
      
c      write(6,200)rpip,rpim,rkp,rkm,rpp,rpb
            
      goto 10

 900  stop
      end

      Subroutine WISER_ALL_SIG(E0,P,THETA_DEG,RAD_LEN,TYPE,SIGMA)

!------------------------------------------------------------------------------
! Calculate pi,K,p  cross section for electron beam on a proton target
! IntegrateQs over function WISER_FIT using integration routine QUADMO
! E0         is electron beam energy, OR max of Brem spectra
! P,E       is scattered particle  momentum,energy
! THETA_DEG  is kaon angle in degrees
! RAD_LEN (%)is the radiation length of target, including internal
!                (typically 5% for SLAC, 2% for JLab)
!               = .5 *(target radiation length in %) +5.
!       ***  =100. IF BREMSTRULUNG PHOTON BEAM OF 1 EQUIVIVENT QUANTA ***
! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
! SIGMA      is output cross section in nanobars/GeV-str
!------------------------------------------------------------------------------

      IMPLICIT NONE       
      REAL E0,P,THETA_DEG,RAD_LEN,SIGMA
      INTEGER TYPE
      COMMON/WISER_ALL/ E,P_COM,COST,P_T,TYPE_COM,PARTICLE ,M_X,U_MAN
      REAL E,P_COM,COST,P_T,M_X,U_MAN
      INTEGER TYPE_COM,PARTICLE
!  Wiser's fit    pi+     pi-    k+     k-     p+      p-   
      REAL A5(6)/-5.49,  -5.23, -5.91, -4.45, -6.77,  -6.53/
      REAL A6(6)/-1.73,  -1.82, -1.74, -3.23,  1.90,  -2.45/
      REAL MASS2(3)/.019488, .2437, .8804/
      REAL MASS(3)/.1396, .4973, .9383/ 
      REAL MP/.9383/,  MP2/.8804/, RADDEG/.0174533/
      REAL  M_L,SIG_E
      REAL*8 E_GAMMA_MIN,WISER_ALL_FIT,QUADMO,E08,EPSILON/.003/
      EXTERNAL WISER_ALL_FIT                        
      INTEGER N,CHARGE
                                            
      P_COM = P
      TYPE_COM = TYPE
      PARTICLE = (TYPE+1)/2       ! 1= pi, 2= K, 3 =P
      CHARGE = TYPE -2*PARTICLE +2  ! 1 for + charge, 2 for - charge
      E08 =E0
                
      E =SQRT(MASS2(PARTICLE) + P**2)

      COST = COS(RADDEG * THETA_DEG)
      P_T = P * SIN(RADDEG * THETA_DEG)
      IF(TYPE.LE.4) THEN  !mesons
       IF(CHARGE.EQ.1) THEN   ! K+ n final state
        M_X = MP
       ELSE   ! K- K+ P final state
        M_X = MP+ MASS(PARTICLE)
       ENDIF
      ELSE  ! baryons 
       IF(CHARGE.EQ.1) THEN   ! pi p  final state
        M_X = MASS(1)  ! pion mass
       ELSE   ! P P-bar  P final state
        M_X = 2.*MP
       ENDIF
      ENDIF
      E_GAMMA_MIN = (M_X**2 -MASS2(PARTICLE ) -MP2+2.*MP*E)/
     >  (2.*(MP -E +P*COST))
!      WRITE(10,'(''E_GAMMA_MIN='',F10.2,''  p_t='',F8.2)')
!     >     E_GAMMA_MIN,P_T
!      E_GAMMA_MIN = MP *(E + MASS(PARTILCE))/(MP -P*(1.-COST))
      
      IF(E_GAMMA_MIN.GT..1) THEN !Kinematically allowed?
       M_L = SQRT(P_T**2 + MASS2(PARTICLE))    

       IF(TYPE.NE.5) THEN  ! everything but proton
        SIG_E = QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N)  *
     >           EXP(A5(TYPE) *M_L) *EXP(A6(TYPE) *P_T**2/E)
       ELSE ! proton

        U_MAN = ABS(MP2 + MASS2(PARTICLE) -2.*MP*E)
        SIG_E = QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N)  *
     >           EXP(A5(TYPE) *M_L) 
       ENDIF
       SIGMA = P**2/E * 1000. * RAD_LEN/100. *SIG_E 
      ELSE   ! Kinematically forbidden
       SIGMA = 0.
      ENDIF

      RETURN
      END

      REAL*8 FUNCTION WISER_ALL_FIT(E_GAMMA)

!---------------------------------------------------------
! Calculates  pi, k, p  cross section for gamma + p -> k
!  It is already divided by E_GAMMA, the bremstrulung spectra
! David Wiser's fit from Thesis, eq. IV-A-2 and Table III.
! Can be called from WISER_SIG using integration routine QUADMO
! E,P are KAON energy and momentum
! P_t is KAON transverse momentum
! P_CM is KAON center of mass momentum
! P_CM_L is KAON center of mass longitudinal momentum
! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
! E_GAMMA is photon energy.
!             Steve Rock 2/21/96
!---------------------------------------------------------
                           
      IMPLICIT NONE       
      COMMON/WISER_ALL/ E,P,COST,P_T,TYPE,PARTICLE,M_X,U_MAN

      REAL  E,P,COST,P_T,M_X,U_MAN
      INTEGER  TYPE  !  1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
      INTEGER PARTICLE   ! 1= pi, 2= K, 3 =P
!  Wiser's fit    pi+     pi-    k+     k-     p+       p- 
      REAL A1(6)/566.,  486.,   368., 18.2,  1.33E5,  1.63E3 / 
      REAL A2(6)/829.,  115.,   1.91, 307.,  5.69E4, -4.30E3 / 
      REAL A3(6)/1.79,  1.77,   1.91, 0.98,  1.41,    1.79 / 
      REAL A4(6)/2.10,  2.18,   1.15, 1.83,   .72,    2.24 /
      REAL A6/1.90/,A7/-.0117/ !proton only
      REAL MASS2(3)/.019488, .2437, .8804/
      REAL MASS(3)/.1396, .4973, .9383/ 
      REAL MP2/.8804/,MP/.9383/, RADDEG/.0174533/
      REAL X_R, S,B_CM, GAM_CM,  P_CM
      REAL P_CM_MAX, P_CM_L
      REAL*8 E_GAMMA
                                            

!Mandlestam variables                                                
      S = MP2 + 2.* E_GAMMA * MP    

!Go to Center of Mass to get X_R
      B_CM = E_GAMMA/(E_GAMMA+MP)
      GAM_CM = 1./SQRT(1.-B_CM**2)
      P_CM_L = -GAM_CM *B_CM *E + 
     >          GAM_CM * P * COST
      P_CM = SQRT(P_CM_L**2 + P_T**2)  


      P_CM_MAX =SQRT (S +(M_X**2-MASS2(PARTICLE))**2/S 
     >    -2.*(M_X**2 +MASS2(PARTICLE)) )/2.
      X_R =  P_CM/P_CM_MAX   
       IF(X_R.GT.1.) THEN  ! Out of kinematic range
        WISER_ALL_FIT = 0.
       ELSEIF(TYPE.NE.5) THEN  ! not the proton
        WISER_ALL_FIT = (A1(TYPE) + A2(TYPE)/SQRT(S)) *
     >   (1. -X_R + A3(TYPE)**2/S)**A4(TYPE)/E_GAMMA  
       ELSE ! special formula for proton
        WISER_ALL_FIT = ( (A1(TYPE) + A2(TYPE)/SQRT(S)) *
     >   (1. -X_R + A3(TYPE)**2/S)**A4(TYPE)          /
     >   (1.+U_MAN)**(A6+A7*S) )/E_GAMMA  
       ENDIF
      
      RETURN
      END

      DOUBLE PRECISION FUNCTION QUADMO(FUNCT,LOWER,UPPER,EPSLON,NLVL)   
      implicit none
      REAL*8 FUNCT,LOWER,UPPER,EPSLON                                   
         INTEGER NLVL                                                   
      INTEGER   LEVEL,MINLVL/3/,MAXLVL/24/,RETRN(50),I                  
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),            
     1   FMRX(50), ESTRX(50), EPSX(50)                                  
      REAL*8  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,       
     1   AREA, ABAREA,   M, COEF, ROMBRG,   EPS                         
         LEVEL = 0                                                      
         NLVL = 0                                                       
         ABAREA = 0.0                                                   
         L = LOWER                                                      
         R = UPPER                                                      
         FL = FUNCT(L)                                                  
         FM = FUNCT(0.5*(L+R))                                          
         FR = FUNCT(R)                                                  
         EST = 0.0                                                      
         EPS = EPSLON                                                   
  100 LEVEL = LEVEL+1                                                   
      M = 0.5*(L+R)                                                     
      COEF = R-L                                                        
      IF(COEF.NE.0) GO TO 150                                           
         ROMBRG = EST                                                   
         GO TO 300                                                      
  150 FML = FUNCT(0.5*(L+M))                                            
      FMR = FUNCT(0.5*(M+R))                                            
      ESTL = (FL+4.0*FML+FM)*COEF                                       
      ESTR = (FM+4.0*FMR+FR)*COEF                                       
      ESTINT = ESTL+ESTR                                                
      AREA=DABS(ESTL)+DABS(ESTR)                                        
      ABAREA=AREA+ABAREA-DABS(EST)                                      
      IF(LEVEL.NE.MAXLVL) GO TO 200                                     
         NLVL = NLVL+1                                                  
         ROMBRG = ESTINT                                                
         GO TO 300                                                      
 200  IF((DABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.                         
     1         (LEVEL.LT.MINLVL))  GO TO 400                            
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0                             
  300    LEVEL = LEVEL-1                                                
         I = RETRN(LEVEL)                                               
         VALINT(LEVEL, I) = ROMBRG                                      
         GO TO (500, 600), I                                            
  400    RETRN(LEVEL) = 1                                               
         MX(LEVEL) = M                                                  
         RX(LEVEL) = R                                                  
         FMX(LEVEL) = FM                                                
         FMRX(LEVEL) = FMR                                              
         FRX(LEVEL) = FR                                                
         ESTRX(LEVEL) = ESTR                                            
         EPSX(LEVEL) = EPS                                              
         EPS = EPS/1.4                                                  
         R = M                                                          
         FR = FM                                                        
         FM = FML                                                       
         EST = ESTL                                                     
         GO TO 100                                                      
  500    RETRN(LEVEL) = 2                                               
         L = MX(LEVEL)                                                  
         R = RX(LEVEL)                                                  
         FL = FMX(LEVEL)                                                
         FM = FMRX(LEVEL)                                               
         FR = FRX(LEVEL)                                                
         EST = ESTRX(LEVEL)                                             
         EPS = EPSX(LEVEL)                                              
         GO TO 100                                                      
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)                          
      IF(LEVEL.GT.1) GO TO 300                                          
      QUADMO = ROMBRG /12.0D0                                           
      RETURN                                                            
      END                                                               



