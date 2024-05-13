         SUBROUTINE H2MODEL(QQ,WW,W1,W2)

********************************************************************************
*
* This subroutine calculates model cross sections for inclusive electron-proton
* scattering in the resonance region. The cross section is returned in 
* nanobarns/sr GeV. This fit is a modified version of Linda Stuart's 8/91 fit. 
* One major difference is that the coefficients represent fit results from a 
* substantially expanded data set, including all inclusive SLAC data in the 
* range 1.1 < W^2 < 5. There are other differences; for a complete discussion, 
* see Keppel's Ph.D. thesis. 2/94 CEK
*
* E        = Incident electron energy in GeV.
* EP       = Scattered electron energy in GeV.
* TH       = Scattered electron angle in degrees.
* SIG1     = Cross section in nb/sr/GeV**2 (DSIG/DOMEGA/DW**2)
* SIG2     = Cross section in nb/sr/GeV (DSIG/DOMEGA/DEP)
* SIG_NRES = Non-resonant contribution to SIG1.
* SIG_RES  = Resonant contribution to SIG1.
* SIG_NR   = Non-resonant contribution to SIG2.
* SIG_R    = Resonant contribution to SIG2.
* SIGroper = Possible Roper contribution to SIG2.
* goroper  = Logical variable, set true if including possible Roper strength
*
* SIG1, SIG2, SIG_NRES, and SIG_RES are 2-parameter arrays where the second
* parameter is the error on the first.
********************************************************************************
        IMPLICIT NONE

        logical goroper 
        logical goodfit 
        INTEGER I
c       REAL*4  E, EP, TH , sige
        REAL*8  SIN2, SIG1(2), SIG2(2), SIG_RES(2), SIG_R(2), 
     >          SIG_RES1(2), SIG_RES2(2), SIG_NRES(2), WW, SIG_NR(2), 
     >          QQ, DIPOLE, COS2, TAU, EPS, K, DW2DEP, DEPCONV, 
     >          PI, AM, ALPHA,W1,W2,NU,CONV,SIGT,SIGL,
     >          RADCON, R_NRES, sigroper(2), x, r, dr,fac,rlog
       
       
        RADCON = 3.141592654/180.
        PI = 3.14159265
        AM = .9382727
        ALPHA = 0.00729735
        CONV = 0.0025767

C        SIN2 = SIN(TH*RADCON/2.0)*SIN(TH*RADCON/2.0)
C        COS2 = 1 - SIN2
C        Q2 = 4.0*E*EP*SIN2
C        W2 = AM*AM + 2.0*AM*(E - EP) - 4.0*E*EP*SIN2

         goroper = .true.

C        write(6,*) q2,w2

        IF(WW.LT.1.15) THEN            ! Below pion threshold
          DO I = 1,2
            SIG1(I) = 0.0
            SIG2(I) = 0.0
            SIG_NRES(I) = 0.0
            SIG_RES(I) = 0.0
            SIG_NR(I) = 0.0
            SIG_R(I) = 0.0
            sigroper(I) = 0.0
          ENDDO
          RETURN
        ENDIF

        K = (WW - AM*AM)/(2.0*AM)
        NU = K + QQ/(2.0*AM)
        TAU = NU*NU/QQ
        x = qq/2./AM/nu
C        EPS = 1.0/(1.0 + 2.0*(1.0 + TAU)*SIN2/COS2)
        DIPOLE = 1.0/(1.0 + QQ/0.71)**2
c        DEPCONV = ALPHA*K*EP/(4.0*PI*PI*Q2*E)*(2.0/(1.0 - EPS))*1000.
c        DW2DEP = 2.0*AM + 4.0*E*SIN2

! H2MOD_FIT returns cross sections in units of microbarns/(dipole FF)**2
        CALL H2MODEL_FIT(QQ,WW,SIG_NRES,SIG_RES1,SIG_RES2,
     >                 sigroper,goroper)
c        call r1990(x,qq,r,dr,goodfit)
c        if (goodfit) then
c         R_NRES = r
c        else 
         R_NRES = 0.25/SQRT(QQ)
c        endif

        SIG_NRES(1) = SIG_NRES(1)*DIPOLE*DIPOLE
C        SIG_NRES(2) = SIG_NRES(2)*DIPOLE*DIPOLE
        SIG_RES(1)  = (SIG_RES1(1) + SIG_RES2(1) + sigroper(1))
     >                 *DIPOLE*DIPOLE
C        SIG_RES(2)  = SQRT(SIG_RES1(2)**2 + SIG_RES2(2)**2 + 
C     >                 sigroper(2)**2)*
C     >                 DIPOLE*DIPOLE
C        SIG1(1) = SIG_NRES(1)*(1.0 + EPS*R_NRES) + SIG_RES(1)
C        SIG1(2) = SQRT( (SIG_NRES(2)*(1.0 + EPS*R_NRES))**2 +
C     >                   SIG_RES(2)**2)
C        SIG2(1) = SIG1(1)*DW2DEP
C        SIG2(2) = SIG1(2)*DW2DEP
C        sig_nr(1) = sig_nres(1)*dw2dep
C        sig_nr(2) = sig_nres(2)*dw2dep
C        sig_r(1) = sig_res(1)*dw2dep
C        sig_r(2) = sig_res(2)*dw2dep
        
C        sige = sig2(1)

         SIGT = SIG_NRES(1)+SIG_RES(1)
c         sigl = r*sigt  
        SIGL = R_NRES*SIG_NRES(1)
         W1 = K*CONV*SIGT/(4.0*ALPHA*PI*PI)
         W2 = K*CONV*(SIGT + SIGL)/(4.0*ALPHA*PI*PI)/(1.0 + TAU)

        RETURN
        END
