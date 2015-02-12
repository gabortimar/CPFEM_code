!***********************************************************************!
!																		!
!        SUBROUTINE ELAST( jgp, el_comp)								!
!-----------------------------------------------------------------------!
!  Routine to give elastic constants (compliance matrix, current axes)	!
!-----------------------------------------------------------------------!
!    Arguments:    jgp		INTEGER; gauss point number					!
!                  el_comp	REAL(6,6); compliance matrix				!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE ELASTC( jgp, el_comp)

	USE FAS_COM			! FASOLT main variables declaration
	USE FAS_CPL			! Crystal plasticity variables


	IMPLICIT NONE
	INTEGER(4)	jgp, iset, jset, i, j, k, l
	REAL(4)		ph1, phi, ph2, a, aug, mult, term1, term2, term3
	REAL(4), DIMENSION(1:6, 1:6):: el_comp
	DIMENSION	a(1:3, 1:3), aug(1:3)

!----------------------------------------------------------current orientation matrix
	ph1= sv(1,jgp) + dsv(1,jgp)
	phi= sv(2,jgp) + dsv(2,jgp)
	ph2= sv(3,jgp) + dsv(3,jgp)

	CALL CTMATRX( ph1, phi, ph2, a, aug)
	a= TRANSPOSE(a)   ! crystal to specimen basis

!-----------------------loop on reduced tensor terms, expand and calculate compliance

	DO iset=1,6
		SELECT CASE(iset)
			CASE(1,2,3)
				i=iset
				j=iset
				mult=1.	    ! multiplier unity if no shear involved 
			CASE(4,5,6)
				i=1+MOD(iset-3,3)
				j=1+MOD(i,3)
				mult=2.		! double if shear involved ( tensor to 'engineering')
		END SELECT

		DO jset=iset,6
			SELECT CASE(jset)
				CASE(1,2,3)
					k=jset
					l=jset
				CASE(4,5,6)
					k=1+MOD(jset-3,3)
					l=1+MOD(k,3)
			END SELECT

				term1=   a(i,1)*a(j,1)*a(k,1)*a(l,1)*e(1,mat_code)	&
					&	+a(i,2)*a(j,2)*a(k,2)*a(l,2)*e(2,mat_code)	&
					&	+a(i,3)*a(j,3)*a(k,3)*a(l,3)*e(3,mat_code)

				term2=	 a(i,2)*a(j,2)*a(k,3)*a(l,3)*e(4,mat_code)	&
					&	+a(i,3)*a(j,3)*a(k,2)*a(l,2)*e(4,mat_code)	&
					&	+a(i,3)*a(j,3)*a(k,1)*a(l,1)*e(5,mat_code)	&
					&	+a(i,1)*a(j,1)*a(k,3)*a(l,3)*e(5,mat_code)	&
					&	+a(i,1)*a(j,1)*a(k,2)*a(l,2)*e(6,mat_code)	&
					&	+a(i,2)*a(j,2)*a(k,1)*a(l,1)*e(6,mat_code)

				term3=	 a(i,2)*a(j,3)*a(k,2)*a(l,3)*e(7,mat_code)	&
					&	+a(i,3)*a(j,2)*a(k,3)*a(l,2)*e(7,mat_code)	&
					&	+a(i,3)*a(j,1)*a(k,3)*a(l,1)*e(8,mat_code)	&
					&	+a(i,1)*a(j,3)*a(k,1)*a(l,3)*e(8,mat_code)	&
					&	+a(i,1)*a(j,2)*a(k,1)*a(l,2)*e(9,mat_code)	&
					&	+a(i,2)*a(j,1)*a(k,2)*a(l,1)*e(9,mat_code)	

			el_comp(iset,jset)= mult*(term1 + term2 + 0.5*term3)
		END DO
	END DO

!----------------------------------------------------------complete symmetric matrix
	DO iset=1,5
		DO jset= iset+1,6
			el_comp(jset,iset)= el_comp(iset,jset)
		END DO
	END DO


END SUBROUTINE ELASTC

!***********************************************************************!
!																		!
!        SUBROUTINE PLASTC( dt, am)										!
!-----------------------------------------------------------------------!
!   Routine to give the slip system resistances and the matrix of		!
!   derivatives of (log)slip resistance WRT slip strain increment.		!
!   Note that the hardening rate is calculated at the current value		!
!   of slip resistance i.e. including effect of slip increments given	!
!-----------------------------------------------------------------------!
!    ARGUMENTS:         dt     REAL; Time step							!
!                       am     REAL; Strain rate sensitivity index		!
!                                              matrix dlnS0/dDgamma		!
!-----------------------------------------------------------------------!
!  Model Parameters:													!
!               ampar(1)   m  slip rate sensitivity						!
!				ampar(2)   theta0 (ds/de at zero stress -thetaIV)		!
!				ampar(3)   thetaIV ( stage IV hardening					!
!				ampar(4)   alpha index									!
!				ampar(5)   tauSS(rate=1)								!
!               ampar(6)   n saturation rate sensitivity				!
!																		!
!***********************************************************************!
SUBROUTINE PLASTC( dt, am)

	USE FAS_COM			! FASOLT main variables declaration
	USE FAS_CPL

	IMPLICIT NONE

	REAL(4)		dt

	INTEGER(4)	i, j, iss, jss, ns
	REAL(4)		THETA, sumss, sterm, arg, sst_iso, sst, ah_1, sums_w, gam_sign, tiny

	INTEGER(4), DIMENSION(1:nss(mat_code)):: iactiv
	REAL(4),	DIMENSION(1:nss(mat_code)):: am, grad_h
	REAL(4),	DIMENSION(1:nss(mat_code),1:nss(mat_code)):: gradient

!---------------------statement function for theta
	THETA(i,arg)= ampar(3,i,mat_code) + &
	& ampar(2,i,mat_code)*SIGN(ABS(arg)**ampar(4,i,mat_code),arg)


!----------------------------------------local control variables
    DATA tiny/ 1.e-14 /

    am(1:nss(mat_code))=ampar(1,1:nss(mat_code),mat_code)

!------------------------------------flag list of active sysytems
	ns= nss(mat_code)

	DO iss= 1, ns
		iactiv(iss)= 0
	END DO

	DO i= 1, nsa
		iss= isa(i)
		iactiv(iss)= i
	END DO

!---------------------------------------accumulated slip increment
	sumss= 0.
	DO iss= 1, ns
		sumss= sumss+ ABS(gam(iss))
	END DO

!-------------------------------------------loop on ALL systems
	DO iss=1, ns
!-------------------------------------rate sensitive target stress
		sst= ampar(5,iss,mat_code)*(sumss/dt)**ampar(6,iss,mat_code)

!-----------------------------------------------current theta
		if(sst .ne. 0.) then
			arg= 1. - ((s0(iss))/sst)
		else
			arg = 0.
		end if

		IF (arg .LT. 0.) arg=0.

		ah_1= THETA(iss,arg)
		grad_h(iss)= ah_1
	END DO

!-------------------------------------------------full hardening matrix
	DO i= 1, ns
		gradient(1:ns,i)= smatrx(1:ns,i,mat_code)*grad_h(i)
	END DO

!--------------------------------------------------slip resistances 
	DO iss= 1, ns
		s1(iss)= s0(iss)
		DO jss=1,ns
			s1(iss)= s1(iss)+ gradient(iss,jss)*ABS(gam(jss))
		END DO
	END DO

!----------------------------actives only for log. hardening matrix
	DO i=1, nsa
		iss=isa(i)
		DO j= 1, nsa
			jss=isa(j)

			IF (ABS(gam(jss)) .GT. tiny) THEN
				gam_sign= gam(jss)/ABS(gam(jss))
				ah(i,j)=  gam_sign* gradient(iss,jss)/s1(iss)
			ELSE
				ah(i,j)= 0.
			ENDIF

		END DO
	END DO

END SUBROUTINE  !PLASTC


!-----------------individual system target(depends on latent matrix)
!		sums_w= 0.
!		DO jss= 1, nss(mat_code)
!			sums_w= sums_w + smatrx(iss,jss,mat_code)*ABS(gam(jss))
!		END DO
!		sst= sst_iso*(1. + sums_w/sumss)

!-----------------------------------------------current theta (RK4)
!		arg= 1. - ((s0(iss)               )/sst) ;IF (arg .LT. 0.) arg=0.
!		ah_1= THETA(iss,arg)

!		arg= 1. - ((s0(iss)+0.5*ah_1*sumss)/sst) ;IF (arg .LT. 0.) arg=0.
!		ah_2= THETA(iss,arg)

!		arg= 1. - ((s0(iss)+0.5*ah_2*sumss)/sst) ;IF (arg .LT. 0.) arg=0.
!		ah_3= THETA(iss,arg)

!		arg= 1. - ((s0(iss)+    ah_3*sumss)/sst) ;IF (arg .LT. 0.) arg=0.
!		ah_4= THETA(iss,arg)

!		grad_h(iss)= (ah_1 + 2.*ah_2 + 2.*ah_3 + ah_4) /6.


