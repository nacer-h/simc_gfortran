	subroutine elastic(e,theta,jmod,csp)
        implicit none

c	q = 4-momentum transfer squared
c	e = incident electron energy

	real*4 mn,nu,mott,x,qq,s2,c2,ef,e,tau,sin2,theta,csp
        real*4 alpha,barn,recoil,w1p,w2p,gep,gmp,t1,t2
        real*4 gen,gmn,QQG,QG,RMUP,gmpcor, gepcor,rhalla,r
        real*4 gepdip,gmpdip
        integer jmod        


	data barn/0.389e-3/ 
        data RMUP/2.792782/
 	alpha = 1/137.03599
        mn = .938272
c	barn: (1 GeV)**-2 = 0.389e-3 barn

        s2=sin(theta/2.)**2
        c2=cos(theta/2.)**2

c	Calculate initial electron energy
        ef = e*mn/(mn+2.*e*s2)
        qq = 4.*e*ef*s2
        nu = e-ef

c        write(6,*) e,ef,theta,qq

 
        tau = qq/4./mn/mn

c	Calculate Mott cross section in barn/(GeV-sr)
	mott = 4*alpha*alpha*ef*ef/qq/qq*c2*barn

c	Calculate recoil factor
	recoil = ef/e

        gep = 1./(1. + qq/0.71)**2 
        gmp = RMUP*gep
         
        gepdip = gep
        gmpdip = gmp

        QQG = qq
        QG = sqrt(qq)
        
        
        if(jmod.EQ.1) then      !!!!  Walker best fit   !!!!

	 QQG=qq

         GEPcor = (1.000-0.23581*QQG+0.51163*QQG**2-0.41838*QQG**3+
     &       0.16334*QQG**4-0.032259*QQG**5+0.0030910*QQG**6-
     &       0.00011312*QQG**7)
         GEP = GEPcor*GEP
         IF(QQG.LT.3.2)GMPcor=(1.000-0.53656*QQG+1.8501*QQG**2-
     &        2.4847*QQG**3+
     &       1.8134*QQG**4-0.79842*QQG**5+0.22068*QQG**6-
     &       0.038443*QQG**7+0.0040854*QQG**8-0.00024114*QQG**9+
     &       0.0000060381*QQG**10)
         IF(QQG.GE.3.2) GMPcor = (1.000+0.027194*QQG-0.0063483*QQG**2+
     &    0.00045490*QQG**3-0.000054447*QQG**4+0.0000032658*QQG**5)
	
 	 GMP = GMPcor*GMP
        
        elseif(jmod.EQ.2) then    !!!!  JRA's fit without Hall A !!!!
  
       
         GMP = RMUP/(1.-0.41031*QG+6.50821*QQG-7.43353*QG*QQG+      !  without Hall A data  ! 
     &         8.03743*QQG**2-2.45134*QG*QQG**2+ 0.40334*QQG**3)
         GEP = 1./(1.+1.24708*QG-4.19625*QQG+15.86942*QG*QQG
     &         -14.50088*QQG**2+7.90098*QG*QQG**2-1.37251*QQG**3)

        elseif(jmod.EQ.3) then    !!!!  JRA's fit with Hall A !!!!

         GMP = RMUP/(1.-0.29940*QG+5.65622*QQG-5.57350*QG*QQG+      !  with Hall A data  !
     &         6.31574*QQG**2-1.77642*QG*QQG**2+0.30100*QQG**3)
         GEP = 1./(1.-0.31453*QG+7.15389*QQG-14.42166*QG*QQG+
     &         23.33935*QQG**2-14.67906*QG*QQG**2+4.03066*QQG**3)
       
       
        elseif(jmod.EQ.4) then    !!!!  Ed Brash's Fit !!!!        

         r = 1.0

         if(QQG.GE.0.04)  r = 1.0 - 0.130*(QQG - .04)
        
         GMP = RMUP/(1. + 0.116*sqrt(qq) + 2.874*qq + 0.241*qq*sqrt(qq) 
     &              + 1.006*qq*qq + 0.345*qq*qq*sqrt(qq))


         GEP = r*GMP/RMUP
        endif

        t1 = tan(theta/2.)
	t2 = t1*t1

	
c	Calculate elastic cross sections for protons 
	w1p = gmp**2 * tau
	w2p = (gep**2 + tau * gmp**2)/(1 + tau)
	
        csp = mott * recoil * (w2p + 2 * w1p * t2)

	return
	end
