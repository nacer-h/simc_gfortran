      SUBROUTINE H2MODEL_FIT(Q2,W2,SIG_NRES,SIG_RES1,SIG_RES2,
     >                    sigroper,goroper)
********************************************************************************
* This is an 24 parameter model fit to NE11 and E133 data and three 
* QSQ points of generated BRASSE data to constrain fits at low QSQ and
* deep inelastic SLAC data from Whitlow to constrain high W2. It has 
* three background terms and three resonances. Each of these is multiplied 
* by a polynomial in Q2. 
*
* 8/91, LMS.
* 7/93. LMS. Modified to include errors.
* SIG_NRES, SIG_RES1, SIG_RES2 are now 2-parameter arrays where the second
* parameter is the error on the first.
********************************************************************************
      IMPLICIT NONE
 
      logical goroper
      INTEGER I, J, KK
      REAL*8 W2, Q2, SIG_NRES(2),SIG_RES1(2),SIG_RES2(2),sigroper(2), 
     >     WR, KLR, KRCM, EPIRCM, PPIRCM, GG, GPI, DEN, 
     >     SIGDEL, W, KCM, K, EPICM, PPICM, WDIF, GAM, 
     >     GR, RERRMAT(25,25), RSIGERR(25),ERRMAT(22,22), SIGERR(22), 
     >     SIGTOTER, ERRCHECK

      REAL*8 PI, ALPHA, AM, MPPI, MDELTA, MPI, GAM1WID, GAM2WID, MASS1,
     >     ROPERWID, MASSROPER, MASS2, DELWID, FIT(30), SQRTWDIF,
     >     XR, MQDEP

       REAL*8 RCOEF(27), COEF(22), RERR1(200),rerr2(125), ERR1(200), 
     >        ERR2(53)
    
      LOGICAL FIRST

      FIRST = .TRUE.

c      pi = 3.14159265
c      alpha = 0.00729735
c      am = .9382727
c      MPPI = 1.07783
c      MDELTA = 1.2340 
c      MPI = 0.13957 
c      GAM1WID = 0.0800 
c      GAM2WID = 0.0900
c      MASS1 = 1.5045 
c      ROPERWID = 0.0500 
c      MASSROPER = 1.4000 
c      MASS2 = 1.6850
c      DELWID = 0.1200
c      XR = 0.1800 
c      MQDEP = 3.40 
      
      pi = 3.14159265
      alpha = 0.00729735
      am = .9382727
      MPPI = 1.07783
      MDELTA = 1.229 
      MPI = 0.13957 
      GAM1WID = 0.080
      GAM2WID = 0.090
      MASS1 = 1.5062 
      ROPERWID = 0.0500 
      MASSROPER = 1.4000 
      MASS2 = 1.6810
      DELWID = 0.120
      XR = 0.1800 
      MQDEP = 3.40 

      DATA RCOEF/
     >   5.2800E+02,  -1.0908E+03,   7.0766E+02,   1.5483E+01, 
     >   4.2450E-01,   8.0152E-01,  -1.9295E+02,   1.0063E+03, 
     >  -6.0730E+02,  -3.7576E+00,   2.8199E+01,   1.8902E+01, 
     >   1.6150E+03,   6.8792E+02,  -1.0338E+03,   2.3285E-01, 
     >   4.6273E-01,   1.7844E-01,   -1.5416E+02, -1.4891E+02,
     >   2.4102E+02,   2.5823E+00,    7.1004E+00, -8.9771E+00,
     >   1.3744E+00,  -1.2085E+00,    1.1218E-01/

      DATA COEF/
     >   4.4050E+02,  -7.9948E+02,   4.8586E+02,   1.5798E+01, 
     >   1.4231E-01,   3.3515E-01,  -2.9657E+02,   1.4930E+03, 
     >  -1.0537E+03,  -3.7598E+00,   2.8633E+01,   1.8381E+01, 
     >   1.6806E+03,   3.2944E+02,  -6.7968E+02,   2.3508E-01, 
     >  -1.6942E+02,  -8.2335E+01,   1.8264E+02,   2.9542E+00, 
     >   5.5004E+00,  -7.7472E+00/
 
c      DATA RERR1/
c     >  2.6120E+02,-9.4211E+02, 4.0844E+03, 7.4994E+02,-3.5317E+03,
c     >  3.1703E+03, 2.1238E-01,-6.1568E-01, 4.1479E-01, 1.9720E-02,
c     >  1.2891E-01,-4.1615E+00, 4.7246E+00, 2.8090E-03, 6.1657E-02,
c     >  1.3120E+00,-9.4379E+00, 9.0902E+00,-1.3636E-03, 2.8054E-02,
c     >  9.6123E-02,-1.1465E+03, 3.9099E+03,-3.0097E+03,-1.0604E+01,
c     > -1.2214E+00,-8.3549E-01, 1.3696E+04, 3.9085E+03,-1.5369E+04,
c     >  1.2745E+04, 2.9942E+01, 7.7268E+00, 1.0113E+01,-4.3868E+04,
c     >  1.5709E+05,-3.0207E+03, 1.2809E+04,-1.1075E+04,-2.0442E+01,
c     > -7.5843E+00,-1.0773E+01, 3.2597E+04,-1.2457E+05, 1.0283E+05,
c     > -1.6960E-01, 5.9410E-01,-4.5486E-01,-1.0715E-02,-2.6512E-03,
c     >  1.0153E-03, 6.4074E+00,-1.9189E+01, 1.3612E+01, 9.2486E-03,
c     >  2.7904E-01, 6.3576E+00,-7.8552E+00, 1.5302E-02,-1.1506E-01,
c     > -4.7552E-02,-1.0171E+01,-1.5884E+00, 1.6223E+01,-1.1379E-04,
c     >  4.9212E-01,-2.4354E+00, 1.7921E+01,-1.7223E+01, 4.0778E-03,
c     > -4.5558E-02,-1.8539E-01, 7.9930E+00,-7.1588E+01, 7.1512E+01,
c     > -2.1529E-03, 1.8337E-01, 7.7590E-01, 7.3007E+02,-2.5219E+03,
c     >  1.9547E+03, 6.1102E+00, 1.2970E+00,-1.3084E+00,-9.4932E+03,
c     >  3.1542E+04,-2.3894E+04,-5.9583E+00, 8.1005E-02, 3.6885E-01,
c     >  9.5708E+03,-2.4911E+03, 9.4342E+03,-7.7120E+03,-1.8608E+01,
c     > -1.1065E+00, 6.5015E+00, 3.1755E+04,-1.1529E+05, 9.1964E+04,
c     >  1.8347E+01,-2.5899E+00, 7.1169E-01,-3.2268E+04, 1.1891E+05,
c     >  1.9339E+03,-7.7737E+03, 6.6128E+03, 1.3392E+01,-7.3587E-02,
c     > -4.9353E+00,-2.4223E+04, 9.2768E+04,-7.6712E+04,-1.3210E+01,
c     >  1.2513E+00,-4.5156E+00, 2.4541E+04,-9.5131E+04, 7.8848E+04,
c     >  1.9952E-02,-7.1332E-02, 5.5522E-02, 9.8804E-04, 2.3682E-04,
c     > -7.9762E-05,-6.3638E-01, 1.9492E+00,-1.4036E+00,-9.9312E-04,
c     > -7.8634E-05, 8.2617E-05, 6.8002E-01,-2.1138E+00, 1.5308E+00,
c     >  1.3008E-04,-1.0389E+02, 3.5942E+02,-2.7883E+02,-6.0671E-01,
c     > -1.3016E-01, 1.4621E-01, 1.2841E+03,-4.3361E+03, 3.3132E+03,
c     >  7.0701E-01, 1.2805E-01, 1.3355E-01,-1.4645E+03, 4.9522E+03,
c     > -3.7686E+03,-1.0047E-01, 2.7406E+02, 3.5483E+02,-1.3433E+03,
c     >  1.0978E+03, 1.9033E+00, 5.3726E-02,-8.1621E-01,-4.3612E+03,
c     >  1.6110E+04,-1.2957E+04,-2.2247E+00,-2.1299E-01,-5.8178E-01,
c     >  4.9755E+03,-1.8393E+04, 1.4724E+04, 3.1774E-01,-9.2555E+02,
c     >  3.4086E+03,-2.7508E+02, 1.1025E+03,-9.3653E+02,-1.4100E+00,
c     >  7.3163E-02, 6.6492E-01, 3.3590E+03,-1.3073E+04, 1.0893E+04,
c     >  1.6311E+00, 2.4826E-01, 8.3308E-01,-3.7999E+03, 1.4772E+04,
c     > -1.2252E+04,-2.3255E-01, 7.0167E+02,-2.7162E+03, 2.2434E+03,
c     >  3.0688E+00,-1.0328E+01, 7.8828E+00, 3.6601E-03, 1.3367E-03,
c     > -2.9672E-03,-3.2441E+01, 1.0979E+02,-8.3795E+01,-6.6345E-03/

c      DATA  rerr2/
c     > 3.7074E-02,
c     >-5.7300E-02, 1.5212E-02, 4.5952E-04,
c     > 1.1568E-04,-2.9315E-04,-4.6018E-01,
c     > 9.3624E-01,-4.5908E-01,-6.2914E-05,
c     > 1.1699E-03, 2.0141E-03, 6.9968E-02,
c     >-1.9348E-01, 1.2176E-01, 5.4214E-07,
c     > 1.3845E-04, 2.5311E-03,-2.5396E-03,
c     >-1.2757E-04, 2.4971E-04,-1.2737E-04,
c     > 7.2023E-03,-4.1404E-03, 4.6704E-04,
c     > -4.6388E-03,-5.2545E-03, 4.0159E+01,-1.3481E+02, 1.0186E+02,
c     >  1.1796E-03,-9.1088E+00, 3.0200E+01,-2.2552E+01, 4.3562E-01,
c     > -1.0404E+01, 3.8414E+01,-3.0978E+01,-1.4730E-02, 4.6327E-03,
c     >  1.9716E-02, 1.1236E+02,-4.1952E+02, 3.3862E+02, 2.4150E-02,
c     >  1.1098E-02, 2.0122E-02,-1.3812E+02, 5.1058E+02,-4.0773E+02,
c     > -4.1791E-03, 3.0702E+01,-1.1132E+02, 8.7622E+01,-1.4199E+00,
c     >  5.0230E+00, 8.0171E+00,-3.1384E+01, 2.6350E+01, 1.3147E-02,
c     > -6.1508E-03,-1.6808E-02,-8.7538E+01, 3.4530E+02,-2.8922E+02,
c     > -1.9581E-02,-1.0895E-02,-2.4705E-02, 1.0611E+02,-4.1369E+02,
c     >  3.4296E+02, 3.2847E-03,-2.3191E+01, 8.8502E+01,-7.2288E+01,
c     >  1.0469E+00,-3.8873E+00, 3.1142E+00,
c     > 1.1348E+00,-1.7657E+00, 4.7686E-01,
c     > 1.6653E-02, 4.3488E-04,-7.5168E-03,
c     >-1.6696E+01, 3.4692E+01,-1.7470E+01,
c     >-4.9697E-03, 4.4232E-02, 5.7617E-02,
c     > 5.7800E+00,-1.3886E+01, 7.9819E+00,
c     > 3.4744E-04,-5.4411E-01, 1.2683E+00,
c     >-7.0771E-01, 1.1282E-02,-2.4800E-02,
c     > 1.2909E-02, 1.5171E-01,-6.0417E-01,
c     > 7.7405E-01,-5.8981E-02,-5.8502E-03,
c     > 8.8611E-04, 5.8326E-03, 6.5418E+00,
c     >-1.2978E+01, 6.1069E+00, 1.2462E-03,
c     >-1.8442E-02,-2.7954E-02,-1.8335E+00,
c     > 4.3674E+00,-2.4393E+00,-6.2354E-05,
c     > 1.4746E-01,-3.4127E-01, 1.8285E-01,
c     >-3.0479E-03, 6.8138E-03,-3.4673E-03,
c     >-7.5270E-02, 4.0914E-02/


c      DATA ERR1/
c     >  3.7797E+02,-1.2732E+03, 4.8470E+03, 9.7589E+02,-3.9592E+03,
c     >  3.3447E+03, 1.9629E-01,-4.2402E-01, 1.9757E-01, 3.0613E-02,
c     > -4.0257E-01,-2.0922E+00, 3.0126E+00, 3.8385E-03, 7.3553E-02,
c     >  1.4084E+00,-8.4718E+00, 7.8586E+00,-1.6484E-03, 2.2185E-02,
c     >  7.4896E-02,-1.5627E+03, 5.0106E+03,-3.7125E+03,-1.1701E+01,
c     > -6.9186E-01,-1.4263E+00, 1.5792E+04, 5.0288E+03,-1.7793E+04,
c     >  1.3974E+04, 3.1643E+01, 5.0040E+00, 9.9958E+00,-4.8540E+04,
c     >  1.6247E+05,-3.7498E+03, 1.4066E+04,-1.1461E+04,-2.0806E+01,
c     > -5.0428E+00,-9.7813E+00, 3.5056E+04,-1.2382E+05, 9.7850E+04,
c     > -2.0038E-01, 5.9769E-01,-4.0397E-01,-1.5776E-02,-3.7509E-03,
c     >  5.7496E-04, 7.2218E+00,-2.0335E+01, 1.3722E+01, 1.2562E-02,
c     >  1.4708E+00, 1.8510E+00,-4.1856E+00, 1.9572E-02,-1.3469E-01,
c     > -3.7791E-02,-1.5215E+01, 1.8843E+01,-9.9384E-01, 5.4133E-04,
c     >  5.6775E-01,-2.4158E+00, 1.5245E+01,-1.4180E+01, 5.3668E-03,
c     > -3.5419E-02,-1.4360E-01, 7.8707E+00,-5.7677E+01, 5.5406E+01,
c     > -7.5727E-04, 1.4127E-01, 5.8964E-01, 1.0277E+03,-3.3407E+03,
c     >  2.4943E+03, 6.1372E+00, 2.0731E+00,-1.0628E-01,-1.1445E+04,
c     >  3.6033E+04,-2.6376E+04,-6.4849E+00,-1.5437E+00,-3.1093E+00,
c     >  1.1966E+04,-3.3062E+03, 1.1473E+04,-8.9323E+03,-1.7658E+01,
c     > -3.0298E+00, 2.4862E+00, 3.6140E+04,-1.2237E+05, 9.3797E+04,
c     >  1.8377E+01, 2.4649E-01, 9.5713E+00,-3.7362E+04, 1.2613E+05,
c     >  2.4733E+03,-8.9836E+03, 7.2301E+03, 1.2133E+01, 1.0120E+00,
c     > -2.0972E+00,-2.6581E+04, 9.4364E+04,-7.4804E+04,-1.2397E+01,
c     >  5.8276E-01,-9.1893E+00, 2.7145E+04,-9.6250E+04, 7.6086E+04,
c     >  2.4070E-02,-7.3772E-02, 5.1165E-02, 1.4597E-03, 3.3977E-04,
c     > -2.6275E-05,-7.2542E-01, 2.0676E+00,-1.4052E+00,-1.3577E-03,
c     > -1.4477E-04,-8.5451E-05, 7.4811E-01,-2.1217E+00, 1.4288E+00,
c     >  1.7439E-04,-1.6022E+02, 5.2231E+02,-3.9172E+02,-4.1771E-01,
c     > -2.3133E-01,-1.9119E-02, 1.6931E+03,-5.4146E+03, 4.0099E+03,
c     >  6.5228E-01, 4.5766E-01, 6.7254E-01,-2.0266E+03, 6.3551E+03,
c     > -4.6404E+03,-9.4689E-02, 4.2768E+02, 5.1531E+02,-1.7829E+03,
c     >  1.3890E+03, 1.1798E+00, 3.1335E-01,-2.5902E-01,-5.3955E+03,
c     >  1.8502E+04,-1.4311E+04,-1.8045E+00,-9.6753E-01,-2.0260E+00,
c     >  6.3626E+03,-2.1445E+04, 1.6387E+04, 2.6350E-01,-1.3456E+03,
c     >  4.5055E+03,-3.8598E+02, 1.3911E+03,-1.1170E+03,-7.9328E-01,
c     > -7.6840E-02, 2.5967E-01, 4.0005E+03,-1.4347E+04, 1.1455E+04,
c     >  1.1795E+00, 6.2629E-01, 1.6961E+00,-4.6485E+03, 1.6399E+04,
c     > -1.2954E+04,-1.7187E-01, 9.8638E+02,-3.4363E+03, 2.7002E+03,
c     >  6.0266E+00,-1.9528E+01, 1.4686E+01,-1.7956E-02, 3.3364E-03,
c     >  1.2080E-03,-5.5018E+01, 1.7933E+02,-1.3517E+02, 7.9955E-03/


c       DATA ERR2/
c     > -2.1546E-02,-2.3493E-02, 7.4315E+01,-2.3518E+02, 1.7398E+02,
c     > -6.4429E-04,-1.9950E+01, 6.3147E+01,-4.6881E+01, 1.2816E+00,
c     > -1.9366E+01, 6.5755E+01,-5.0971E+01, 5.7005E-02, 3.3439E-04,
c     >  5.5786E-03, 1.7715E+02,-6.1369E+02, 4.7999E+02,-2.9558E-02,
c     >  5.5461E-02, 7.1075E-02,-2.3560E+02, 7.9154E+02,-6.0792E+02,
c     >  2.7242E-03, 6.3265E+01,-2.0981E+02, 1.6050E+02,-4.0749E+00,
c     >  1.3388E+01, 1.4562E+01,-5.1058E+01, 4.0557E+01,-4.3474E-02,
c     > -4.4868E-03,-6.3041E-03,-1.3274E+02, 4.7814E+02,-3.8441E+02,
c     >  2.5677E-02,-3.8538E-02,-5.8204E-02, 1.7424E+02,-6.0799E+02,
c     >  4.8014E+02,-2.6425E-03,-4.6992E+01, 1.6058E+02,-1.2570E+02,
c     >  3.0554E+00,-1.0258E+01, 7.9929E+00/


! Kinematic variables.
      IF(FIRST) THEN
        FIRST = .FALSE.
        KLR = (MDELTA*MDELTA - AM*AM)/(2.0*AM)
        KRCM = KLR*AM/SQRT(AM*AM + 2.0*KLR*AM)
        EPIRCM = 0.5*(MDELTA*MDELTA + MPI*MPI - AM*AM)/MDELTA
!Define error matrix:
c        KK = 0
c        if (goroper) then 
c          DO J = 1,25
c            DO I = 1,J
c              KK = KK + 1
c              if (kK.le.200) RERRMAT(I,J) = RERR1(KK)
c              if (kK.le.325.and.kK.gt.200) RERRMAT(I,J)=RERR2(KK-200)
c            ENDDO
c          ENDDO
c        endif
c       if (.not.goroper) then  
c          DO J = 1,22
c            DO I = 1,J
c              KK = KK + 1
c              if (kK.le.200) ERRMAT(I,J) = ERR1(KK)
c              if (kK.le.253.and.kK.gt.200) ERRMAT(I,J)=ERR2(KK-200)
c            ENDDO
c          ENDDO
c        endif
c        if (goroper) then
c          DO J = 1,25
c              DO I = J+1,25
c                RERRMAT(I,J) = RERRMAT(J,I)
c              ENDDO
c          ENDDO
c        endif
c        if (.not.goroper) then
c          DO J = 1,22
c              DO I = J+1,22
c                ERRMAT(I,J) = ERRMAT(J,I)
c              ENDDO
c          ENDDO
c        endif
      ENDIF

      PPIRCM = SQRT(MAX(0.0,(EPIRCM*EPIRCM - MPI*MPI)))
      W = SQRT(W2) 
      WDIF = MAX(0.0001,W - MPPI)
      K = (W*W - AM*AM)/(2.0*AM)
      EPICM = (W*W + MPI*MPI - AM*AM)/(2.0*W)
      PPICM = SQRT(MAX(0.0,(EPICM*EPICM - MPI*MPI)))
      KCM = K*AM/SQRT(AM*AM + 2.0*K*AM)
      GG = DELWID*(KCM/KRCM)**2*(KRCM*KRCM + XR*XR)/
     >     (KCM*KCM + XR*XR)
      GPI = DELWID*(PPICM/PPIRCM)**3*
     >      (PPIRCM*PPIRCM + XR*XR )/(PPICM*PPICM + XR*XR)
      DEN = (W*W - MDELTA*MDELTA)**2 + (MDELTA*GPI)**2
      SIGDEL = 389.4*2.0*PI*ALPHA*(W/AM)*(KRCM/KCM)*
     >         (Q2/K)*GG*GPI/DELWID/DEN

! Get each of the components of the model. 
      SQRTWDIF = SQRT(WDIF)
      FIT(1) = SQRTWDIF
      FIT(2) = WDIF*SQRTWDIF
      FIT(3) = WDIF*WDIF*SQRTWDIF
      FIT(4) = SIGDEL
      if (goroper) FIT(25) = ROPERWID/((W - MASSROPER)**2 + 
     >         0.25*ROPERWID*ROPERWID)
      FIT(5) = GAM1WID/((W - MASS1)**2 + 0.25*GAM1WID*GAM1WID)
      FIT(6) = GAM2WID/((W - MASS2*(1.0 + Q2*MQDEP/1000.0))**2 + 
     >         0.25*GAM2WID*GAM2WID)
      DO I = 1,6
        FIT(I + 6)  = FIT(I)*Q2
        FIT(I + 12) = FIT(I)*Q2*Q2
      ENDDO
      DO I = 1,3
        FIT(I + 18)  = FIT(I)*Q2*Q2*Q2
        FIT(I + 21)  = FIT(I)*Q2*Q2*Q2*Q2
      ENDDO
      if (goroper) FIT(26)  = FIT(25)/sqrt(Q2)
      if (goroper) FIT(27)  = FIT(25)/Q2

! Find sig_t (in microbarns/gd**2).
      SIG_NRES(1) = 0.0
      SIG_RES1(1) = 0.0
      SIG_RES2(1) = 0.0
      SIG_NRES(2) = 0.0
      SIG_RES1(2) = 0.0
      SIG_RES2(2) = 0.0
      SIGTOTER = 0.0
      SIGroper(1) = 0.0
      SIGroper(2) = 0.0
      if (goroper) then
        DO J = 1,27
c          RSIGERR(J) = 0.0
c          DO I = 1,25
c            RSIGERR(J) = RSIGERR(J) + FIT(J)*FIT(I)*RERRMAT(I,J)
c            SIGTOTER = SIGTOTER + FIT(J)*FIT(I)*RERRMAT(I,J)
c          ENDDO
          IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12.OR.J.EQ.17
     > .OR.J.EQ.18 ) THEN
             SIG_RES2(1) = SIG_RES2(1) + FIT(J)*RCOEF(J)          
c             SIG_RES2(2) = SIG_RES2(2) + RSIGERR(J)
          ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
             SIG_RES1(1) = SIG_RES1(1) + FIT(J)*RCOEF(J)
c             SIG_RES1(2) = SIG_RES1(2) + RSIGERR(J)
          elseIF(j.ge.25.and.j.le.27) then
            SIGroper(1) = SIGroper(1) + FIT(J)*RCOEF(J)          
c            SIGroper(2) = SIGroper(2) + RSIGERR(J)
          ELSE
            SIG_NRES(1) = SIG_NRES(1) + FIT(J)*RCOEF(J)
c            SIG_NRES(2) = SIG_NRES(2) + RSIGERR(J)
          ENDIF
        ENDDO
      endif
      if (.not.goroper) then
        DO J = 1,22
c          SIGERR(J) = 0.0
c          DO I = 1,22
c            SIGERR(J) = SIGERR(J) + FIT(J)*FIT(I)*ERRMAT(I,J)
c            SIGTOTER = SIGTOTER + FIT(J)*FIT(I)*ERRMAT(I,J)
c          ENDDO
          IF(J.EQ.5.OR.J.EQ.6.OR.J.EQ.11.OR.J.EQ.12) THEN
             SIG_RES2(1) = SIG_RES2(1) + FIT(J)*COEF(J)          
c             SIG_RES2(2) = SIG_RES2(2) + SIGERR(J)
          ELSEIF(J.EQ.4.OR.J.EQ.10.OR.J.EQ.16) THEN
            SIG_RES1(1) = SIG_RES1(1) + FIT(J)*COEF(J)
c            SIG_RES1(2) = SIG_RES1(2) + SIGERR(J)
          ELSE
            SIG_NRES(1) = SIG_NRES(1) + FIT(J)*COEF(J)
c            SIG_NRES(2) = SIG_NRES(2) + SIGERR(J)
          ENDIF
        ENDDO
      endif

! ERRCHECK should agree with SIGTOTER.
C      ERRCHECK = SQRT(ABS(SIG_RES2(2) + SIG_RES1(2) + SIG_NRES(2)))
C      SIGTOTER = SQRT(SIGTOTER)
c      SIG_RES2(2) = SQRT(ABS(SIG_RES2(2)))
c      SIG_RES1(2) = SQRT(ABS(SIG_RES1(2)))
c      SIG_NRES(2) = SQRT(ABS(SIG_NRES(2)))

      RETURN
      END
