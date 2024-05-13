      PROGRAM RATES
      IMPLICIT none

      integer*8 i,j,n,nbins
      real*4 epass,eppass,thpass,thpass1,thpass2,sigpass
      real*8 sigmac,sigelast,sigma1,sigma2,sigma3
      real*8 sigmael,pi,mp,mp2,m,m2,e,ep,theta,de
      real*8 nu,q2,w2,epsilon,sin2,sigma,theta1,delmax,delmin
      real*8 bcur,lum,accpt,k1temp,rateinel,rate(50)
      real*8 trateinel,tinc,thetamin,epin,sigmaave
      real*8 radcon,tan2,dens,tlength,sigmatot,ratetot
      real*8 w2min,w2max,dw2,q2c,w2c,epelast,rateel
      logical elast,sos,shms,fileok
      character*80 infile
      character*3 sos_flag

      write(6,'('' input file: '',$)')
      read(5,'(a)') infile

      open(30,file=infile,access='sequential',status='old',err=20)
      fileok = .true.
      sos=.false.
      shms=.false.
 20   continue
      if (.not.fileok) then
         write(6,*)' Error opening input file'
         stop
      endif

      open(31,file='rates.out',access='sequential',status='new')
 
      read(30,'(a)') sos_flag
      if (sos_flag.eq.'SOS') then
         sos = .true.
         shms= .false.
         write(6,*)' Rates will be for SOS '
      elseif (sos_flag.eq.'HMS') then
         sos = .false.
         shms= .false.
         write(6,*)'   Rates will be for HMS '
         write(31,*)'   Rates will be for HMS '
      elseif (sos_flag.eq.'SHM') then
         sos = .false.
         shms= .true.
         write(6,*)'   Rates will be for SHMS '
         write(31,*)'   Rates will be for SHMS '
      else 
         write(6,*)' Unrecognized sos_flag '
         stop
      endif

c      write(6,*)sos,shms
c gh: ad-hoc 1.6x adjustment to solid angle is to make rate estimates closer to SIMC
c gh: 18.10.11 - widen to full momentum acceptance for rates, not fiducial
      if (sos.eqv..true.) then 
         delmax = 0.25
         delmin = 0.25
         accpt = 9.0e-3*1.6
c         write(6,*)' SOS '
      elseif (shms.eqv..true.) then    ! SHMS
         delmax = 0.20
         delmin = 0.20
          accpt = 4.0e-3*1.6 
c         write(6,*)' SHMS '
      elseif (sos.eqv.shms) then  ! HMS
         delmax = 0.10
         delmin = 0.10
         accpt = 6.5e-3*1.6
c         write(6,*)' HMS '
      endif
      nbins = int(delmax*100.+delmin*100.)
c      write(6,*)' del',sos,shms,delmin,delmax,nbins
      
      pi = 3.141592654
      radcon = pi/180.
      mp = .9382727
      mp2 = mp*mp
      dens = 0.0723
      tlength = 10.0
      rateinel = 0.
      elast = .false.
      
      write(31,30)
 30   format(t7,"E",t18,"E'",t25,"theta",t35,"I",t42,"q2",t48,"w2"
     *     ,t53,"epsilon",t62,"sig_el",t71,"el_rate",t81,"sig_inel",
     *     t91,"inel_rate",t101,"tot_rate",/
     *     t5,"(GeV)",t16,"(GeV)",t25,"(deg)",t34,"(uA)",t62,
     *     "(nb/sr)",t72,"(kHz)",t80,"(nb/sr/GeV)",t93,"(kHz)",
     *     t102,"(kHz)")

c      write(6,*) "Enter E, E', theta, current"
 40   read(30,*,end=950) e,epin,theta,bcur            

      epass = e
      eppass = epin
      thpass = theta*radcon
      lum = 3758.6*bcur*dens*tlength     
      call model(epass,eppass,thpass,11,sigpass)
      sigmac = sigpass
      tinc = 28.*2./20./1000./radcon
      thetamin = theta - 10.5*tinc
      sigmatot = 0.
      de = 0.01*epin  
      sin2 = sin(radcon*theta/2.)*sin(radcon*theta/2.)
      q2c = 4.*e*epin*sin2
      nu = e-epin
      w2c = mp2+2.*nu*mp-q2c
      
CCC   Find out whether run contains elastics and calculate   CCC
CCC   elastic scattering rates if true.                      CCC
      
      epelast = e*mp/(mp+2.*e*sin2)   
      if(((1.+delmax)*epin).GE.epelast) then 
         elast = .true.
         thpass1 = thpass - 28./1000. - 28.*2./20/1000.
         sigmael = 0.
         do j = 1,20
            thpass1 = thpass1 + 28.*2./20/1000.
            call elastic(epass,thpass1,2,sigpass)
            sigmael = sigmael + sigpass
         enddo
         sigmael = 1.0e+9*sigmael/20.
 1       rateel = accpt*lum*sigmael
      else
         rateel = 0.
      endif
      
c      write(6,*) epin,epelast,elast,sigmael,rateel
      
CCC   Calculate inelastic rates                              CCC
      
      ep = epin - (100.*delmin-0.5)*de
      rateinel = 0.
      do i=1,nbins              !!!  Loop over Energy bins and sum rates !!!
         ep = ep + de
         nu = e-ep
         sigmaave = 0.
         do j=1,20              !!!  Loop over angle bins and average rates  !!!
            theta = thetamin + j*tinc
            sin2 = sin(radcon*theta/2.)*sin(radcon*theta/2.)
            tan2 = sin2/(1.-sin2)
            q2 = 4.*e*ep*sin2
            w2 = mp2+2.*nu*mp-q2
            epsilon = 1./(1.+2.*(nu*nu/q2+1.)*tan2)
            epass = e
            eppass = ep
            thpass = theta*radcon
            if(w2.GT.1.18) then
               call model(epass,eppass,thpass,11,sigpass)
            else
               sigpass = 0.
            endif
            sigmaave = sigmaave + sigpass*1.d0
         enddo
         sigmaave = sigmaave/20.
         rate(i) = accpt*lum*de*sigmaave
         rateinel = rateinel + rate(i)
         
c     write(6,*) i,rate(i),rateinel,rateel
         
      enddo
      ratetot = rateinel + rateel
      
c      write(6,*) rateinel,rateel,ratetot
      
      ratetot = ratetot/1000.   !!!  kHz

      write(31,920) e,epin,theta,bcur,q2c,w2c,epsilon,
     *     sigmael,rateel/1000.,sigmaave,rateinel/1000.,
     *     ratetot
      
 920  format(3f10.4,f7.1,3f7.3,5e10.3)

      goto 40

 950  close(30)
      close(31)
      
      end
