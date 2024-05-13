       SUBROUTINE H2MODEL_IOANA(QSQ,WSQ,W1,W2)
********************************************************************************
*
* This subroutine calculates model cross sections for H2 in resonance region.
* Cross section is returned in nanobarns. This model is valid in the QSQ
* range 0.75 - 10.0 (GeV/c)**2 and the WSQ range from pion threshold to
* 3.0 GeV.
*
* QSQ      = 4-MOMENTUM TRANSFER SQUARED (GEV/C)**2
* WSQ      = Missing mass squared (GeV)**2
* W1,W2    = Inelastic structure functions. (1/GeV)
*
* 8/91, LMS.
* 2/93, LMS: Modified to return W1 and W2 structure functions instead of
*       cross sections. For original version of H2MODEL go to
*       $2$DIA3:[OFLINE.INELAS].
*****************************************************************************
        IMPLICIT NONE

        INCLUDE 'imodel.cmn'

        REAL*8  WSQ, QSQ
        REAL*8  W1, W2,  R_NRES, SIG_RES(2), SIG_RES1(2),
     >          SIG_RES2(2), SIG_NRES(2), SIGT, SIGL, DIPOLE, K, NU,
     >          TAU, PI, AM, ALPHA,
     >          CONV,sig_roper(2)
	integer i
	logical goroper
c	logical first
*
        INTEGER IHD 
*
        COMMON/tmptgt/ IHD

        goroper = .false.
c	first = .true.
 
*       
*
* N.I. ...
*
        DATA PI/3.14159265/, AM/.938259/, ALPHA/0.00729735/,
     >          CONV/0.0025767/
*
*
* if first time around than do yourself a favour and read the parameters
*
	if (first) then
           write(6,*)'Reading parameters for Ioana''s model'
           first = .not.first
           open(unit=48,file
     >          ='parm_h2.txt'
     >          ,status='old')
           do i=1,36
              read(48,*)j
              read(48,*)xvalh(i)
              read(48,*)tmptmp
              read(48,*)tmptmp
              read(48,*)tmptmp
           enddo	
           close(48)
        endif
*
! Check that kinematic range is valid.
        W1 = 0.0
        W2 = 0.0
        IF(WSQ.LT.1.15.OR.WSQ.GT.4.3.OR.
     >    QSQ.GT.10.0) THEN
           do i=1,2	
              SIG_NRES(i) = 0.0
              SIG_RES1(i) = 0.0
              SIG_RES2(i) = 0.0
              SIG_ROPER(i)= 0.0 
           enddo
*          WRITE(6,'(''  H2MODEL_IOANA called outside of kinematic range'')')
          RETURN
        ENDIF

!  returns transverse cross sections in units of
! microbarns/(dipole FF)**2
        CALL i_h2_model(QSQ,WSQ,SIG_NRES,SIG_RES1,SIG_RES2,
     &                 SIG_ROPER,goroper,xvalh)

        NU = (WSQ + QSQ - AM*AM)/(2.0*AM)
        TAU = NU*NU/QSQ
        K = (WSQ - AM*AM)/(2.0*AM)
        DIPOLE = 1.0/(1.0 + QSQ/0.71)**2
        R_NRES = 0.25/SQRT(QSQ)           ! Corresponds to R used in fits.
        SIG_NRES(1) = SIG_NRES(1)*DIPOLE*DIPOLE
        SIG_RES(1)  = (SIG_RES1(1) + SIG_RES2(1)+sig_roper(1))*
     &                 DIPOLE*DIPOLE
        SIGT = SIG_NRES(1) + SIG_RES(1)
        SIGL = R_NRES*SIG_NRES(1)
! The factor CONV converts from GeV*microbarns to 1/GeV
        W1 = K*CONV*SIGT/(4.0*ALPHA*PI*PI)
        W2 = K*CONV*(SIGT + SIGL)/(4.0*ALPHA*PI*PI)/(1.0 + TAU)
        RETURN
        END
