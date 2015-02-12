!***********************************************************************!
!																		!
!                       Fasolt3   -cpstrs								!
!                        Pete Bate 2000									!
!-----------------------------------------------------------------------!
!            Crystal plasticity constituitive routine for 3-d FEM		!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!

!***********************************************************************!
!																		!
!     SUBROUTINE CPSTRS( jgp, el_comp, de, r, s, t, ep, chi)			!
!-----------------------------------------------------------------------!
!     Routine to give the elasto-plastic stress state for crystal		!
!     plasticity by two levels of iteration, given a reasonable			!
!     estimate of stress state to start with. The plastic strain		!
!     increments and the hardening modulus are returned for use in		!
!     the e-p Jacobian contribution. State variable increments are		!
!     created, on a temporary basis, in this routine and downstream		!
!-----------------------------------------------------------------------!
!    Arguments:    jgp		INTEGER; Gauss point number					!
!                  el_comp	REAL ARRAY(6,6); elastic compliances		!
!                  de		REAL ARRAY(6); Strain increments			!
!                  r		REAL ARRAY(3); R.B. rotation about x,y,z	!
!                  s		REAL ARRAY(6); Final stresses				!
!                  t		REAL ARRAY(6), Initial stress (current axes)!
!                  ep		REAL ARRAY(6); Plastic strain increments	!
!                  chi		REAL ARRAY(6,6); Elasto-plastic modulus		!
!																		!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE CPSTRS( jgp, el_comp, del, rl, s, t, epl, chi)

	USE FAS_COM			! FASOLT main variables declaration
	USE FAS_CPL			! Crystal plasticity variables

	IMPLICIT NONE

	INTEGER(4)	jgp	! index of current IP

	LOGICAL(4)	l_iter
	INTEGER(4)	i, j, k, kk, iter, ifail, maxcpsit, iop
	REAL(4)		dt, ph1, phi, ph2, capfac, gamcap, scale, erepmd, &
				&	FJ2VOM, sum, dp1, dph, dp2, tolercp, alpha_cps

    	REAL(4), DIMENSION( 1:3,1:3):: a
	REAL(4), DIMENSION(1:6):: del, epl, s, t, erreps
	REAL(4), DIMENSION(1:3):: aug, rl
	REAL(4), DIMENSION(1:6,1:6):: el_comp, chi, djac
	REAL(4), DIMENSION(1:maxss):: am

!---------------------------------------Common for sum of iteration error
	REAL(4)		sumcpiter
	COMMON/ CPITER / sumcpiter
    
!--Local control data: capfac is max. slip/eff. strain inc, alpha and max. iterations.
	DATA capfac/  6.0/  !4

	alpha_cps= alpha
	maxcpsit=  maxit
	tolercp=   toler

!-----------------Pull down time step, orientation and slip resistances
    dt=tiem(inc)
    ph1= sv(1,jgp) + dsv(1,jgp)
    phi= sv(2,jgp) + dsv(2,jgp)
    ph2= sv(3,jgp) + dsv(3,jgp)
    	DO i= 1, nss(mat_code)
		s0(i)= sv(i+3,jgp)
		s1(i)= s0(i) + dsv(i+3,jgp)
	END DO

    	gamcap= capfac*FJ2VOM(.FALSE., del)	!	Maximum slip allowed
	CALL CTMATRX(ph1,phi,ph2,a,aug)		!	Orientation matrix 
    	CALL BTMATRX(a)						!	Slip matrix

!----------------------------------------------Iterate for stress state
	iter=0
	l_iter= .TRUE.
	DO WHILE( l_iter )
		iter=iter+1
!-----------------------------------------------------------Slip stresses 
        	tau(1:nss(mat_code))= MATMUL(TRANSPOSE(b(1:6,1:nss(mat_code))),s)
			
!----------------------------------------Slip increment predictor with cap
		CALL SLPINI(jgp,dt,am)
		CALL SCALE_IT(am,s,gamcap,scale)

!-------------------Slip increment refinement, resistances and visco-plastic matrix
		CALL SLPFTA( dt)
		CALL SCALE_IT( am, s, gamcap, scale)

!------------------------------Current guess at plastic strain inc.	
		epl= MATMUL(b(1:6,1:nss(mat_code)),gam(1:nss(mat_code)))
	
		erreps= del - epl -MATMUL( el_comp, ( s - t) )  ! Form the strain inc. error

!-----------------------------Factor ah matrix, form Jacobian
		CALL FACTOR(nsa,ah(1:nsa,1:nsa),ifail)
		IF ( ifail .EQ. 1 ) CALL FAFERR(910)
		DO i= 1, 6
			DO j= 1, 6
				DO k= 1, nsa
					kk= isa(k)
					wks(k)= b(j,kk)
				END DO
				CALL SOLVER(nsa,ah(1:nsa,1:nsa),wks(1:nsa))
				sum=0.
				DO k= 1, nsa
					kk= isa(k)
					sum= sum + b(i,kk)*wks(k) 
				END DO
				djac(i,j)= sum
			END DO
		END DO

		djac= djac + el_comp	!	Add elastic compliances

!--------------------------------------------------------Correction to stress
		CALL FACTOR(6,djac,ifail)
		IF (ifail.EQ.1) CALL FAFERR(910)
		CALL SOLVER(6,djac,erreps)

		erepmd= FJ2VOM(.TRUE.,erreps)/FJ2VOM(.TRUE.,s) ! Stress error measure

!------------------------------------Damping to prevent 'fly-off'
		scale= 1./( 1. + alpha_cps * erepmd )

		DO i=1,6
			erreps(i)= scale*erreps(i)
		END DO
!-------------------------------------Correct the stress
		DO i=1,6
			s(i)= s(i) + erreps(i)
		END DO

		l_iter = (iter.LT.maxcpsit) .AND. (erepmd.GT.tolercp) ! Test for repeat

	END DO
	! Write slip rates to srate variable:
	do iop=1,nss(mat_code)
		srate(iop, jgp) = ABS(gam(iop) / dt)
	end do
	! (this will be overwritten several times because all this is in the outermost NR loop,
	! but this is the only way...)

	sumcpiter= sumcpiter + erepmd				! CP iteration error for report

	CALL ECHROT(a,aug,rl,dp1,dph,dp2)  ! Estimate of orientation change

!-----------------------Calculate elasto-plastic modulus (invert) 
	DO i=1,6
		erreps(1:6)=0.
		erreps(i) = 1.
		CALL SOLVER(6,djac,erreps)
		DO j=1,6
			chi(j,i)= erreps(j)
		END DO
	END DO

!		 Put current orientation change and slip resistance
!					 change terms into state variable list
    dsv(1,jgp)=dp1
    dsv(2,jgp)=dph
    dsv(3,jgp)=dp2
	DO i=1,nss(mat_code)
		dsv(i+3,jgp)=s1(i)-s0(i)
	END DO

END SUBROUTINE CPSTRS

!***********************************************************************!
!																		!
!     SUBROUTINE SCALE_IT( am, gamcap, scale)							!
!-----------------------------------------------------------------------!
!	     Routine to scale slip stress and strains increments			!
!-----------------------------------------------------------------------!
!    Arguments:		am		REAL;		Slip rate sensitivity			!
!					gamcap	REAL;			Cap on magnitude of slips	!
!					scale	REAL;			Scaling for slips			!
!																		!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE SCALE_IT(  am, s, gamcap, scale)

	USE FAS_COM			! FASOLT main variables declaration
	USE FAS_CPL			! Crystal plasticity variables

	IMPLICIT NONE
	REAL(4)		amax, gamcap, scale
	INTEGER(4)	i, ii, k
	REAL(4)		gammax, smult, SPOWER
	REAL(4), DIMENSION(1:maxss):: am
	REAL(4), DIMENSION(1:6):: s

!------------------get biggest slip increment and its rate sensitivity
	amax= am(1)
	gammax=ABS(gam(isa(1)))

	DO i= 2, nsa
		ii= isa(i)
		IF (ABS(gam(ii)) .GT. gammax) THEN
			gammax= ABS(gam(ii))
			amax= am(ii)
		END IF
	END DO

!------------------calculate approximate scaling for strain and stress
	IF (gammax .GT. gamcap) THEN
		scale= gamcap/gammax
		smult= SPOWER(0.9*scale, afactor*amax)

! This is a fiddle for small m values only and may never happen!
				IF (smult.EQ.0.) smult=0.004/amax

		DO i= 1, 6
			s(i)= smult*s(i)
		END DO
		DO i= 1, nss(mat_code)
			gam(i)= scale*gam(i)
			tau(i)= smult*tau(i)
		END DO
	ELSE
		scale=1.
	ENDIF
       
END SUBROUTINE SCALE_IT
