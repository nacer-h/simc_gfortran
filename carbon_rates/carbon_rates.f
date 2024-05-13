      program xem2_rates_times_signuc2

      implicit none

      real aa,zz,nn,esep
      common /nuc/aa,zz,nn,esep

      real ebs,ths,eps
      real sigma

      real*8 eprime
      real*8 atmp,blo,bhi,q2lo,q2hi,pelo,pehi,pe,deltape,deleit
      real*8 sigtot,sigit
      real*8 Pmin,Pmax
      real*8 xlo,xhi,eplo,ephi
      real*8 xmin,xmax
      real*8 eb,thdeg,deltap,Pcent
      real*8 q2min,q2max
      real*8 thstart,pstart
      real*8 A,Z
      real*8 mp
      real*8 thr,cs,sn,tn
      real*8 Q2,nu,WSQ,x
      real*8 F1,F2,Rc,W1,W2
      real*8 sigmott,sig
      real*8 Na ! Avagoddro's number
      real*8 qe
      real*8 dOmega,luminosity,lrad
      real*8 thick,cur,Ntarg,Ne,rate
      real*8 t40k,t10k
      real*8 ps,tmax,ttmp,xtmax,targ_time
      integer nset
      integer i,j,nxbin,k

      parameter(mp=0.938272)
      parameter(dOmega=6.0d-3) ! 4 msr for SHMS, 6 (?) for HMS
      parameter(Na=6.02214129d23)
      parameter(qe=1.602176565d-19)



      write(6,*) 'Beam energy (GeV):'
      read(5,*) eb

      write(6,*) 'Spectrometer momentum (GeV): '
      read(5,*) pcent

      write(6,*) 'Spectrometer angle (deg.): '
      read(5,*) thdeg

      write(6,*) 'Beam current (uA): '
      read(5,*) cur

      cur=cur/1.0d6

** Inputs
c 10 cm hydrogen
c      A=1.0
c      Z=1.0
c      thick = 10.0*0.0723
c      Ntarg=Na*thick/1.00794    !atoms/cm2
c      aa=1.0
c      zz=1.0
c      nn=0.0
c      esep=0.0
c      cur=60.0d-6

C 10 cm deuterium
c      A=2.0
c      Z=1.0
c      thick = 10.0*0.167
c      Ntarg=Na*thick/2.014    !atoms/cm2
c      aa=A
c      zz=z
c      nn=aa-zz
c      esep=2.1/1000.0
c      cur=60.0d-6


C carbon target, 1.5%
      A=12.0
      Z=6.0
      thick=0.1749
      Ntarg=Na*thick/12.0107 !atoms/cm2
      aa=A
      zz=z
      nn=aa-zz
      esep=17.27/1000.0
c      cur=60.0d-6      


 
c      cur=15.0d-6
      Ne=cur/qe ! electrons/sec

      luminosity=Ne*Ntarg


      thr=thdeg*3.141592654/180.0
      cs = cos(thr/2.)
      sn = sin(thr/2.)
      tn = tan(thr/2.)

      write(6,*) 'Beam energy', eb, ' GeV'
      write(6,*) 'Scattering angle',thdeg, ' degrees'
      write(6,*) 'Pcent', Pcent
        
      pelo=0.9*pcent
      pehi=1.1*pcent

      deleit=(pehi-pelo)/100.0
      q2 = 4.*eb*pcent*sn**2
      nu = eb - pcent
      x = q2/2./mp/nu
      wsq=-q2+mp**2+2*mp*nu

      write(6,*) 'Other kinematics: '
      write(6,*) 'Q2: ', q2
      write(6,*) 'W2: ', wsq
      write(6,*) 'x: ', x
      
      sigtot=0.0
      do j=1,100                !sum over pe for each x-bin
         eprime = pelo + real(j-1.)*deleit
         nu = eb - eprime
         q2 = 4.*eb*eprime*sn**2
         x = q2/2./mp/nu
         wsq=-q2+mp**2+2*mp*nu

         ebs=eb
         ths=thr
         eps=eprime
         sigit=sigma(ebs,ths,eps)*1000.0 !nb/sr/GeV
C stuff for using F1F2IN09
c     sigint = sigint + sigtot*delta*1000
c            if (wsq.gt.0) then 
c               call F1F2IN09(Z, A, Q2, Wsq, F1, F2, Rc)
c            else
c               F1=0.0
c               F2=0.0
c            endif
C     Convert F1,F2 to W1,W2
c            W1 = F1/mp
c            W2 = F2/nu
c nb/GeV/sr
c            sigmott=(0.197320/(2.0*137.0388*eb*sn**2))**2*cs**2*1.0d7
c            sigit = sigmott*(W2+2.0*W1*tn**2)
C**********************************
         sigit = sigit*deleit
         sigtot = sigtot + sigit
      enddo                     ! loop over pe
c     C need to get rates now.
      ps=1.0
      rate=sigtot*dOmega*luminosity*1d-33

      write(6,*) 'Integrated cross section: ',sigtot
      write(6,*) 'Total electron rate (Hz): ',rate


c

c

c
c
c            
c         write(4,33) x,thdeg,ep,wsq,Q2,sig,rate,t40k,t10k
c      enddo

 22	format(12(3x,a10))
 33	format(7(3x,f10.6),3x,E10.5,3x,4(3x,f8.2))
 44	format(10(3x,a7))
c 55    format(f10.4,1x,f10.2,1x,f9.3,2x,f10.3,2x,f10.2,3x,e10.4,3x,e10.4,
c     >       3x,e10.4,3x,e10.4)
      end



      function sigma(ei,the,ep)
      real aa,zz,nn,esep
      common /nuc/aa,zz,nn,esep
      w = ei - ep
      sin2 = (sin(the/2))**2
      q2 = 4.*ei*ep*sin2
      x = q2/2./.938/w
      recoil = ei/( 1. + 2.*ei*sin2/.94/aa)
      if (ep.ge.recoil) then
      	sdeep = 0.
      	sqe = 0.
      	go to 10
      endif
      call ftosig(ei,the,w,zz,nn,esep,sigma,y)
      if (sigma.lt.1e-20) sigma = 0.
      sqe = sigma
      izz = zz
      inn = nn
      sdeep = 0
      if (x.ge.1.) go to 10
      call struct(ei,ep,sin2,sdeep,smott,2,izz,inn,aa,vw2)
C...sdeep is in pb/sr/GeV; want nb/sr/MeV
      sdeep = sdeep*1e-6
  10  sigma = sqe + sdeep
*      print14,ep,sdeep,sqe,sigma
      write(33,*) vW2
      write(33,14)ep,sdeep,sqe,sigma
14    format(f7.3,3(E17.3))
      return
      end

      subroutine ftosig(e1,th,w,zz,nn,eb,sigma,y)

      real e1,th,w,zz,nn,eb,sigma,y
c
c
c	Subroutine to compute sigma (nb/sr/MeV)
c		from f(y) in 1/GeV
c	th is theta in radians
c	e1,w in GeV
c	y returned in GeV/c
c	eb is binding energy in GeV
c
c
      e2 =e1 - w
      tr = th
      t2 = tr/2.
      s2 = sin(t2)**2
      q2 = 4.*e1*e2*s2
      tt2 = tan(t2)**2
      a = zz+nn
      am = .938*a
      tau = q2/4./.938**2
      call hoehler(q2,gep,gmp,gen,gmn)
      ge2 = zz*gep**2+nn*gen**2
      gm2 = zz*gmp**2+nn*gmn**2
      w20 = (gm2*tau+ge2)/(1.+tau)
      w10 = tau*gm2
C     print*,'all OK'
      y = ycalc(eb,w,q2,a,ierr,dwdy)
      q3 = sqrt(q2+w**2)
      pf = (1.-exp(-a/8.))*.22+.04
      pp2 = .4*pf**2
      g2 = sqrt(y**2+pp2+((a-1.)*.938)**2)
      rec = 1 + 2.*e1*s2/.938
      sigm = cos(t2)/(2.*137.*e1*s2)
      sigm = (.01973*sigm)**2/rec
      qrat = q2/q3**2
      fact2 = ((am - g2 - w*y/q3)**2 + (1. - w**2/q3**2)*pp2/2.)/.938**2
      w2 = fact2*w20/dwdy
      w1 = (w10 + w20*pp2/2./.938**2)/dwdy
      fact = sigm*(w2 + 2.*w1*tt2)
      f = real(fy(y,a))
      sigma = 1.e6*f*fact
      sigep = 1.e6*sigm*(w20 + 2.*w10*tt2)
      if (ierr.eq.1) sigma = 0.
      return
      end

      function ycalc(esep,eloss,q2,aa,ierr,dwdy)
      real*8 m,m2,a,b,c,z,v,w,beta,q3
      ierr = 0
C     print*,'OK here too'
      m = .938
      w = eloss
      m2 = m**2
      beta = aa*m + w - esep
      pf = (1.-exp(-aa/8.))*.22+.04
      pp2 = .4*pf**2
      q3 = dsqrt(w**2 + q2)
      z = q3**2 + m2 - (aa - 1)**2*m2 - beta**2
      a = 4.*(q3**2 - beta**2)
      b = 4.*q3*z
      c = z**2 - 4.*(pp2 + m2*(aa - 1.)**2)*beta**2
      v = b**2 - 4.*a*c
      if (v.lt.0.) then
      v = 0.
      ierr = 1
      endif
      y = - b - dsqrt(v)
      y = y/2./a
      if (y.gt..1) ierr = 2
      f1 = 4*q3*w*y*y - 2*beta*y*y + 2*z*w*y - 2*aa*q3*y - aa*m*z
      f1 = f1 - 2*beta*pp2 - 2*beta*(aa - 1)**2*m*m
      f2 = 2*beta**2*y - 2*q3*q3*y - z*q3
C     print*,f1,f2
      if (f1.eq.0.) then
      ierr = 1
      dwdy = 0
      else
      dwdy = f2/f1
      endif
      ycalc = y
      return
      end

      real function fy(y,a)
        real y
      if (a.gt.4) then
        pf=.26
      else
      pf=(1.-exp(-a/8.))*.22+.04
      end if
      y0=pf/1.8
      gy=(y/y0)**2
      yy=11.513*y
      fy=.290*exp(-gy/2.)/y0+2./(exp(yy)+exp(-yy))
      return
      end

      subroutine hoehler(t,gep,gmp,gen,gmn)
c
c	Computes nucleon form factors according to
c	Hoehler et al. NP B114, 505(1976)
c	t=Q2>0 in GeV/c **2
c
      f1r=(.955+.09/(1.+t/.355)**2)/(1.+t/.536)/2.
      f2r=(5.335+.962/(1.+t/.268))/(1.+t/.603)/2.
      tw=.783**2
      fw=.67/(tw+t)
      f1p=f1r+fw+(-.39)/(.9216+t)+(-.54)/(2.756+t)
      f2p=f2r+(.04)/(tw+t)+(-1.88)/(1.30+t)+.24/(10.18+t)
      tphi=1.02**2
      f1v=f1r+.05/(1.464+t)+(-.52)/(6.+t)+.28/(8.703+t)
      f2v=f2r+(-1.99)/(1.464+t)+.2/(6.+t)+.19/(8.703+t)
      f1s=.71/(tw+t)+(-.64)/(tphi+t)+(-.13)/(3.24+t)
      f2s=(-.11)/(tw+t)+(.13)/(tphi+t)+(-.02)/(3.24+t)
      f1n=f1s-f1v
      f2n=f2s-f2v
      gep=f1p-t*f2p/4./.88035
      gmp=f1p+f2p
      gen=f1n-t*f2n/4./.8828
      gmn=f1n+f2n
      return
      end

      subroutine STRUCT(E0,EP,SINSQ,O1,O2,IFLAG,IZTARG,intarg,
     >     AATARG,vw2)

C...**************************************************************
C...THIS ARE THE COEFFICIENTS NEEDED TO CORRECT FOR NEUTRON EXCESS
C...AND THE EMC EFFECT. THEY WERE OBTAINED BY FITTING THE RATIOS OF
C...CROSS SECTIONS [CS(A>2)/CS(D2)] OBTAINED IN E139.THE CORRECTING
C...FUNCTION IS A 6-TH ORDER POLYNOMIAL IN X.

      REAL EF(7) /-0.00136693,-0.00510425,-0.0375986,-0.0946004,
     >  -0.122435,-0.0112751,0.406435/

C...**************************************************************

C     COMMON/TGTNUC/ IZTARG, INTARG, AATARG, RLTARG, TGTDEN, ITNAME,
C    >               IZWALL, INWALL, AAWALL, RLWALL, WALDEN, IWNAME


      GOTO (1,2) IFLAG

C...***************************************************************
C...THIS IS THE FIT USED FOR THE INELASTIC IN E139 WITHOUT CORRECTION
C...FOR THE EMC EFFECT.
C...THE PROCESS IS AS FOLLOWS: (1) OBTAIN THE INELASTIC CROSS
C...SECTION FOR D2 USING BODEK FIT, (2) CALCULATE THE CROSS
C...SECTION FOR A NUCLEUS WITH (A,Z',N'), WHERE Z'= NUMBER OF
C...PROTONS, N'= NUMBER OF NEUTRONS. HERE A=Z'+N' AND Z'=N'.
C...(3) CORRECT FOR NEUTON EXCESS TO OBTAIN THE CROSS SECTION FOR
C...A NUCLEUS WITH (A,Z,N).

C...GET THE CROSS SECTION FOR D2 USING BODEK FIT.

1     CALL XSECHD(E0,EP,SINSQ,1,2,VW2,W2,W1,XSECTN,XMOTT,X)

C...HERE WE CALCULATE THE CROSS SECTION FOR (A,Z',N')

      XSECMD= XSECTN * AATARG/2.

C...HERE WE CORRECT FOR NEUTRON EXCESS AND OBTAIN THE FINAL CROSS
C...SECTION.
C THE APPROXIMATION 1-.8X IS USED FOR SIGN/SIGP (SEE BODEK PHYS REV
C D, 1979) WHICH MAY HAVE PROBLEMS FOR X>.8 BUT IS PRETTY GOOD
C BELOW THAT (UNCERTAINTY IN SIGN/SIGP IS ~+/-.04 FOR ALL X).

      RNP=1.-.800*ABS(X)
      ZA = FLOAT(IZTARG)/AATARG
      CNEUT = 0.5*(1.+RNP)/(ZA+(1.-ZA)*RNP)
      O1 = XSECMD/CNEUT
      O2 = XMOTT
      RETURN
C...***************************************************************

C...***************************************************************
C...THIS IS THE FIT USED FOR THE INELASTIC IN E139 WITH CORRECTION
C...FOR THE EMC EFFECT.
C...THE PROCESS IS AS FOLLOWS: (1) OBTAIN THE INELASTIC CROSS
C...SECTION FOR D2 USING BODEK FIT, (2) CALCULATE THE CROSS
C...SECTION FOR A NUCLEUS WITH (A,Z',N'), WHERE Z'=NUMBER OF
C...PROTONS, N'=NUMBER OF NEUTRONS. HERE A=Z'+N' AND Z'=N'.
C...(3) CORRECT FOR NEUTRON EXCESS AND DIFFERENCE BETWEEN
C...D2 AND A NUCLEUS WITH A>2 BY USING FIT TO THE RATIOS OF
C...THEIR CROSS SECTION AS OBTAINED IN E139.

C...GET THE CROSS SECTION FOR D2 USING BODEK FIT.

2     if ( aatarg.gt.1.) then
          CALL XSECHD(E0,EP,SINSQ,1,2,VW2,W2,W1,XSECTN,XMOTT,X)
C         ...HERE WE CALCULATE THE CROSS SECTION FOR (A,Z',N')
          XSECMD = XSECTN * AATARG/2.
      else
          CALL XSECHD(E0,EP,SINSQ,1,1,VW2,W2,W1,XSECTN,XMOTT,X)
          XSECMD = XSECTN
      endif


C...HERE WE CORRECT FOR NEUTRON EXCESS AND THE EMC EFFECT.

      SUMEF=EF(1)
      DO 100 J=2,7
      ZZ=J-1.
100   SUMEF=SUMEF+EF(J)*X**(ZZ)
      RNTD2=1.+SUMEF*ALOG(AATARG)
      IF(aaTARG.EQ.1) RNTD2=1.0
      O1 = XSECMD * RNTD2
      O2 = XMOTT
      RETURN
C...**************************************************************

      END


      SUBROUTINE XSECHD(E0,EP,SINSQ,ITYP,IHD,VW2,W2,W1,XSECTN,CSMOTT,
     >X)
C...
C...THIS FITS ARE FOR HYDROGEN AND DEUTERIUM.
C...USES ATWOOD'S RESONANCE FITS AND CALCULATE VW2 AS:
C...VW2=B*F2 WHERE B IS BACKGROUND WITH RESONANCES AND F2 IS
C...A UNIVERSAL FUNCTION TO THE DEEP INELASTIC
C...
C...PARAMETER DESCRIPTION:
C...E0 IS THE ELECTRON INCIDENT ENERGY.                  (INPUT)
C...EP IS THE SCATTERED ELECTRON ENERGY.                 (INPUT)
C...SINSQ IS (SIN(THETA/2))**2.THETA IS THE SCATTERING
C...ANGLE                                                (INPUT)
C...ITYP CHOOSES WHICH FIT IS GOING TO BE USED.          (INPUT)
C...IHD= 1 FOR HYDROGEN,IHD= 2 FOR DEUTERIUM.
C...ONLY VALID IF THE FIT IS AVAILABLE FOR H OR D2       (INPUT)
C...VW2,W2,W1 AND XSECTN ARE THE FITTED STRUCTURE
C...FUNCTIONS AND THE FITTED CROSS SECTION RESPECTEVILY. (OUTPUT)
C...XMOTT IS THE MOTT CROSS SECTION                      (OUTPUT)
C...X IS THE PARTON MOMENTUM FRACTION, X = QSQ/2MV       (OUTPUT)
C...
      REAL*8 C(24),CF(11),CD(24),CFD(11),F2,B,QSQ,WW
      DATA PMSQ/.880324/,TPM/1.876512/
      DATA R/.18/,ALPHAX/137.0388/
C...
C...::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C...CROSS SECTION FIT FOR HYDROGEN OR DEUTERIUM
C...USES ATWOOD'S RESONANCE FIT WITH COEFFICIENTS FROM FITTING
C...E87,E49A,E49B.COEFFICIENTS SUPLIED BY ARIE BODEK FOR E139
C...FINAL HYDROGEN COEFFS. FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE
C...CHISQ=3995,NO. OF POINTS=2533,FREE PARAMS=28,CHISQ/D.F.=1.59
C...THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)
C...
C...C(24)=HYDROGEN COEFFICIENTS FOR B (BACKGROUND AND RESONANCE
C... TERMS)
C...
      DATA   C(1) / 0.10741163D 01/,  C(2) / 0.75531124D 00/,
     *       C(3) / 0.33506491D 01/,  C(4) / 0.17447015D 01/,
     *       C(5) / 0.35102405D 01/,  C(6) / 0.10400040D 01/,
     *       C(7) / 0.12299128D 01/,  C(8) / 0.10625394D 00/,
     *       C(9) / 0.48132786D 00/,  C(10)/ 0.15101467D 01/,
     *       C(11)/ 0.81661975D-01/,  C(12)/ 0.65587179D 00/,
     *       C(13)/ 0.17176216D 01/,  C(14)/ 0.12551987D 00/,
     *       C(15)/ 0.74733793D 00/,  C(16)/ 0.19538129D 01/,
     *       C(17)/ 0.19891522D 00/,  C(18)/-0.17498537D 00/,
     *       C(19)/ 0.96701919D-02/,  C(20)/-0.35256748D-01/,
     *       C(21)/ 0.35185207D 01/,  C(22)/-0.59993696D 00/,
     *       C(23)/ 0.47615828D 01/,  C(24)/ 0.41167589D 00/
C...
C...CF(11)=HYDROGEN COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION)
C...OMEGAW FIT
C...
      DATA CF(1) / 0.25615498D 00/,  CF(2) / 0.21784826D 01/,
     *     CF(3) / 0.89783738D 00/,  CF(4) /-0.67162450D 01/,
     *     CF(5) / 0.37557472D 01/,  CF(6) / 0.16421119D 01/,
     *     CF(7) / 0.37635747D 00/,  CF(8) / 0.93825625D 00/,
     *     CF(9) / 0.10000000D 01/,  CF(10)/ 0.0           /,
     *     CF(11)/ 0.50000000D 00/
C...
C...
C... FINAL DEUTERIUM COEFFS FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE
C... CHISQ=4456,NO. OF POINTS 2303,FREE PERAMS=26,CHISQ/D.F.=1.96
C... THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)
C...
C... CD(24)=DEUTERIUM COEFFICIENTS FOR B (BACKGROUND AND RESONAN
C... TERMS)
C...
      DATA  CD(1) / 0.10521935D 01/, CD(2) / 0.76111537D 00/,
     *      CD(3) / 0.41469897D 01/, CD(4) / 0.14218146D 01/,
     *      CD(5) / 0.37119053D 01/, CD(6) / 0.74847487D 00/,
     *      CD(7) / 0.12399742D 01/, CD(8) / 0.12114898D 00/,
     *      CD(9) / 0.11497852D-01/, CD(10)/ 0.14772317D 01/,
     *      CD(11)/ 0.69579815D-02/, CD(12)/ 0.12662466D 00/,
     *      CD(13)/ 0.15233427D 01/, CD(14)/ 0.84094736D-01/,
     *      CD(15)/ 0.74733793D 00/, CD(16)/ 0.19538129D 01/,
     *      CD(17)/ 0.19891522D 00/, CD(18)/-0.24480414D 00/,
     *      CD(19)/ 0.14502846D-01/, CD(20)/-0.35256748D-01/,
     *      CD(21)/ 0.35185207D 01/, CD(22)/-0.21261862D 00/,
     *      CD(23)/ 0.69690531D 01/, CD(24)/ 0.40314293D 00/
C...
C...CFD(11) ARE DEUTERIUM COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION)
C....OMEGAW FIT
C...
      DATA CFD(1) / 0.47708776D 00/, CFD(2) / 0.21601918D 01/,
     *     CFD(3) / 0.36273894D 01/, CFD(4) /-0.10470367D 02/,
     *     CFD(5) / 0.49271691D 01/, CFD(6) / 0.15120763D 01/,
     *     CFD(7) / 0.35114723D 00/, CFD(8) / 0.93825625D 00/,
     *     CFD(9) / 0.10000000D 01/, CFD(10)/ 0.0           /,
     *     CFD(11)/ 0.50000000D 00/
C...::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C...
C...COMPUTE SOME KINEMATIC QUANTITIES
C...
      COSSQ=1.0-SINSQ
C...
C...CHECK THAT DIVISION BY ZERO DOES NOT HAPPEN
C...
      IF((E0.EQ.0.).OR.(EP.EQ.0.).OR.(SINSQ.EQ.0.).OR.(
     * COSSQ.EQ.0.)) GOTO 4
C...
      TANSQ=SINSQ/COSSQ
      Q2=4.0*E0*EP*SINSQ
      V=E0-EP
      X= Q2/(TPM*V)
      VSQ=V*V
      W=SQRT(PMSQ+TPM*V-Q2)
      CSMOTT=(19732./(2.0*ALPHAX*E0*SINSQ))**2*COSSQ
      OMEGAP=TPM*V/Q2+PMSQ/Q2
C...
C...OVERCOME RISK OF UNDERFLOW IN THE EXPONENTIATION
C...
      OMEGAP=AMIN1(20.0,OMEGAP)
C...
      SP=1.0-EXP(-7.7*(OMEGAP-1.0))
      WW=W
      QSQ=Q2

C...CHOOSE THE FIT
      GOTO (1) ITYP

C...THIS IS ATWOOD'S TYPE OF FIT WITH COEFFICIENTS OBTAINED FROM
C...FITTING E87,E49A AND E49B.
C...CHOOSE IF HYDROGEN OR DEUTERIUM FIT

1     IF(IHD.EQ.2) GOTO 2

C...GET UNIVERSAL AND RESONANCE FIT FOR HYDROGEN
C...
      UNIV=F2(WW,QSQ,CF)
      BRES=B(WW,QSQ,C)
      GOTO 3
C...
C...GET UNIVERSAL AND RESONANCE FIT FOR DEUTERIUM
C...
   2  UNIV=F2(WW,QSQ,CFD)/SP
      BRES=B(WW,QSQ,CD)
C...
C...COMPUTE VW2,W2,W1,AND THE CROSS SECTION
C...
   3  VW2=UNIV*BRES
      W2=VW2/V
      W1=(1.0+VSQ/QSQ)/(V*(1.0+R))*VW2
      XSECTN=(W2+2.0*TANSQ*W1)*CSMOTT
      RETURN
   4  XSECTN=0.
      RETURN
      END
C...
C...
      DOUBLE PRECISION FUNCTION F2(WM,QSQ,CF)
C...
C... UNIVERSAL FUNCTION FOR ATWOOD'S FIT
C...
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION CF(11)
      DATA PM2/1.876512D0/,PMSQ/.880324D0/
C...
C... OMEGAW FIT...NO PHOTO-PRODUCTION COUPLING
C...
      V=(WM**2+QSQ-PMSQ)/PM2
      OMEGA=2.D0*CF(8)*V/QSQ
      XX=1.D0/OMEGA
      XPX=CF(9)+CF(10)*(XX-CF(11))**2
      OMEGAW=(2.D0*CF(8)*V+CF(6))/(QSQ+CF(7))
      ARG=1.D0-1.D0/OMEGAW
      F2=OMEGAW/OMEGA*ARG**3*(CF(1)+CF(2)*ARG+
     > CF(3)*ARG**2+CF(4)*ARG**3+CF(5)*ARG**4)
      F2=F2*XPX
      RETURN
      END


      DOUBLE PRECISION FUNCTION B(WM,QSQ,C)
C...
C...BACKGROUND AND RESONANCE CONTRIBUTION FOR ATWOOD'S FIT
C...
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(24)
      DIMENSION LSPIN(4)
      DATA LSPIN/1,2,3,2/
      DATA PMSQ/.880324D0/,PM2/1.876512D0/,PM/.938256D0/
      DATA NRES/4/,NBKG/5/
C...
C...KINEMATICS
C...
      WSQ=WM**2
      OMEGA=1.D0+(WSQ-PMSQ)/QSQ
      X=1.D0/OMEGA
      XPX=C(22)+C(23)*(X-C(24))**2
      PIEMSQ=(C(1)-PM)**2
C...
C...COLLECT BACKGROUND TERMS
C...CHECK FOR EXPONENTIAL UNDERFLOWS BEFORE THEY HAPPEN
C...
      B1=DMAX1(0.D0,(WM-C(1)))/(WM-C(1))*C(2)
      EB1=C(3)*(WM-C(1))
      IF(EB1.GT.25.0D0) GO TO 1
      B1=B1*(1.0D0-DEXP(-EB1))
   1  B2=DMAX1(0.D0,(WM-C(4)))/(WM-C(4))*(1.D0-C(2))
      EB2=C(5)*(WSQ-C(4)**2)
      IF(EB2.GT.25.0D0) GOTO 2
      B2=B2*(1.0D0-DEXP(-EB2))
   2  CONTINUE
      BBKG=B1+B2
      BRES=C(2)+B2
C...
C...COLLECT RES. CONTRIBUTION
C...
      RESSUM=0.D0
      DO 30 I=1,NRES
      INDEX=(I-1)*3+1+NBKG
      RAM=C(INDEX)
      IF(I.EQ.1) RAM=C(INDEX)+C(18)*QSQ+C(19)*QSQ**2
      RMA=C(INDEX+1)
      IF(I.EQ.3) RMA=RMA*(1.D0+C(20)/(1.D0+C(21)*QSQ))
      RWD=C(INDEX+2)
      QSTARN=DSQRT(DMAX1(0.D0,((WSQ+PMSQ-PIEMSQ)/(2.D0*WM))**2-PMSQ))
      QSTARO=DSQRT(DMAX1(0.D0,((RMA**2-PMSQ+PIEMSQ)/
     > (2.D0*RMA))**2-PIEMSQ))
      IF(QSTARO.EQ.0.D0) GOTO 40
      TERM=6.08974D0*QSTARN
      TERMO=6.08974D0*QSTARO
      J=2*LSPIN(I)
      K=J+1
      GAMRES=RWD*(TERM/TERMO)**K*(1.D0+TERMO**J)/(1.D0+TERM**J)
      GAMRES=GAMRES/2.D0
      BRWIG=GAMRES/((WM-RMA)**2+GAMRES**2)/3.1415926D0
      RES= RAM*BRWIG/PM2
      GOTO 30
 40   RES=0.D0
 30   RESSUM=RESSUM+RES
C...
C...FORM VW2/F2
C...
      B=BBKG*(1.D0+(1.D0-BBKG)*XPX)+RESSUM*(1.D0-BRES)
      RETURN
      END
