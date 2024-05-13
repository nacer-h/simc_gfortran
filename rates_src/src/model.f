      SUBROUTINE MODEL(epass,eppass,thpass,TARGET,SIG4)

      implicit none
        

**************************************************************************************
*                                                                                    *
*   Calculates and returns model cross section.  H2(D2)model is used for W2 < w2min, *
*   f2glob is used for W2 > w2max, and an interpolation between the two is done      *
*   for w2min < W2 < w2max.  For solid targets f2glob is used above the resonance    *
*   region and a quasielastic model is added for the resonance region.               *
*                                                                                    *
*   Energies are passed in GeV, theta is passed in radians.  Sigma is returned in    *
*   nb/Sr/Gev.                                                                       *  
**************************************************************************************  

      INCLUDE 'flags.cmn'

      real*8 E,EP,TH
      real*8 WW,QQ,x,xtemp,nu,F2DIS,F2RES,amu(20)
      real*8 W1,W2,MP,MP2,mott,frac,q4,sig,emcres,emcdis
      real*8 w1res,w2res,sin2,cos2,tan2,alpha,alpha2,dr,qqtemp
      real*8 w2max,w2min,rres,rdis,sigres,sigdis
      real*4 epass,eppass,thpass,sig4
      integer tar,target, taremc
      character*1 tarf2
      logical goodfit

      usenmc = .false.  !!! only use SLAC f2 param for now  !!!

      sig = 0.

      e = epass
      ep = eppass
      th = thpass

      w2max = 4.3
      w2min = 3.5

CCCCCC            Target to Amu conversion table                    CCCCCCC

      amu(2) = 63                   ! Cu
      amu(3) = 63                   ! Cu
      amu(4) = 12                   ! C12
      amu(5) = 12                   ! C12
      amu(11) = 1               
      amu(15) = 2
      amu(17) = 27                    

CCCCCC        Change from external target #s to internal ones       CCCCCCC

      if(target.EQ.15) then
       tar = 2                         !   Deuterium   !
      elseif(target.EQ.11) then
       tar = 1                         !   Hydrogen    !
      else
       tar = 10
      endif                          !!!   Heavy targets need to be included  !!!

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(tar.EQ.1) tarf2='H'           !   For input into F2glob  !
      if(tar.GE.2) tarf2='D'
      mp = .9382727
      mp2 = mp*mp
      alpha = 1/137.03599
      alpha2 = alpha*alpha  
      sin2 = sin(th/2.)*sin(th/2.)
      cos2 = cos(th/2.)*cos(th/2.)
      tan2 = tan(th/2.)*tan(th/2.)
      qq = 4.*e*ep*sin2
      qqtemp = qq
      if(qqtemp.LT.0.35) qqtemp = 0.35
      q4 = qq*qq
      nu = e-ep
      ww = mp2 + 2.*mp*nu -qq
      x = qq/2./mp/nu


c      write(6,*) "MOD WW:  ",ww

      mott = (4*alpha2*ep*ep*cos2/q4)*0.389e6        !  In nb/GeV/sr  !   

      if(target.EQ.11.OR.target.EQ.15) taremc = 1                    !!! H2 or D2  !!!
      if(target.EQ.2.OR.target.EQ.3) taremc = 7                      !!!  Cu       !!!         
      if(target.EQ.4.OR.target.EQ.5) taremc = 5                      !!!  CARBON   !!!          
      if(target.EQ.17) taremc = 6                                    !!! ALUMINUM  !!!
      if(target.EQ.1) taremc = 9                                     !!!   GOLD    !!!

      if(WW.LE.w2max) then       !!!!  Calculate model for Resonance Region  !!!!                  

       emcres = 1.
       if(tar.EQ.1) then
         call h2model(qq,ww,w1,w2)
         call h2model_liang(qq,ww,w1,w2)
         w1res = w1
         w2res = w2
       endif 

       if(tar.EQ.2) then 
        call d2model_ioana(qq,ww,w1,w2)        
        w1res = w1/2.
        w2res = w2/2. 
       endif

       if(tar.GT.2) then 
        if(x.GT.0.95) x = 0.99           !!!!   Make sure it doesn't blow up  !!!!
        if(usenmc) then                  !!!!   Use Antje's Parameterization of F2   !!!!
         call f2allm(x,qq,f2res)
         f2res=f2res*0.5*(0.976+x*(-1.34+x*(1.319+x*
     &          (-2.133+x*1.533)))+1)
        else
         call f2glob(x,qq,tarf2,9,f2res)    !!! change later to include     QE  !!!
        endif

        call r1990(x,qqtemp,rres,dr,goodfit)

        if(rres.LE.0.) write(6,*) "Warning:  R <= 0."

        w2res = f2res/nu
        w1res = (1.0 + nu*nu/qq)*w2/(1.0 + rres)
        call emcparm(taremc,x,emcres)
        w2res = w2res*emcres
       endif                      
              
       sigres = mott*(w2res + 2.*w1res*tan2)

c       write(6,*) "MODEL:  ",x,qqtemp,sigres

      endif                           !!!!  End Resonance Region  !!!!


      if(WW.GE.w2min) then            !!!!   Calculate Model for DIS   !!!!        

       if(usenmc) then  
        call f2allm(x,qq,f2dis)       !!!!  Use Antje's Parameterization of F2  !!!!
        if(tar.GE.2) f2dis=f2dis*0.5*(0.976+x*(-1.34+x*(1.319+x*
     &          (-2.133+x*1.533)))+1)
       else
        call f2glob(x,qq,tarf2,9,f2dis) !!! change later to correct target !!!
       endif
      
       call emcparm(taremc,x,emcdis)
       f2dis = f2dis*emcdis
       call r1990(x,qqtemp,rdis,dr,goodfit)
       xtemp = x

       w2 = f2dis/nu 
       w1 = (1.0 + nu*nu/qq)*w2/(1.0 + rdis)
       
       if(rdis.EQ.0.) write(6,*) "Warning:  R = 0."
        
       sigdis = mott*(w2 + 2.*w1*tan2)

c       write(6,*) "MODEL2:  ",ww,sigdis

      endif                           !!!!   End DIS   !!!!

      if(WW.GT.w2min.AND.WW.LT.w2max) then
       frac = (ww - w2min)/(w2max - w2min)
       sig = sigdis*frac + sigres*(1.0 - frac)
      elseif(WW.LE.w2min) then
       sig = sigres
      elseif(WW.GE.w2max) then
       sig = sigdis
      endif

      sig4 = sig

      end














