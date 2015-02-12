!***********************************************************************!
!																		!
!                       Fasolt3   -crys_cons							!
!                        Pete Bate 2000									!
!-----------------------------------------------------------------------!
!               Crystal plasticity routines for 3-d FEM					!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!

!***********************************************************************!
!																		!
!      SUBROUTINE SLPFTA( dt)											!
!-----------------------------------------------------------------------!
!   Routine for the slip increments for a given orientation and slip	!
!   stresses, by Newton-Raphson on the log.-log. basis. This version	!
!	tests for the absolute change in slips and returns the augmented	!
!	visco-plastic hardening matrix (X)									!
!-----------------------------------------------------------------------!
!   Arguments:	dt		REAL;			 Time step						!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE SLPFTA( dt)

	USE FAS_COM			! FASOLT main variables declaration
	USE FAS_CPL			! Crystal plasticity variables

	IMPLICIT NONE

	INTEGER(4)	jgp
	REAL(4)		dt

	LOGICAL(4)	l_iter
	INTEGER(4)	i, j, ii, jj, iter, ifail
	REAL(4)		facmax, tiny, conv, gamn, taufg
	REAL(4)		tauerr, errmax, scale, SPOWER

	REAL(4), DIMENSION(1:nss(mat_code)):: am
	REAL(4), DIMENSION(1:nsa):: errtau

    DATA facmax, tiny/  5., 1.E-12/			!	Other data

!---------------------------------------------Tolerance for slip stress convergence
    conv= toler*MAXVAL(ampar(1,1:nss(mat_code),mat_code))	

!------------------------------------------------------Iteration for slip increments
    iter=0
	l_iter= .TRUE.

	DO WHILE( l_iter)
		iter= iter + 1
		CALL PLASTC( dt, am)	! Hardening update

!------------------------------------Complete log. slip ah matrix and factor
		DO i= 1, nsa
			ii=isa(i)
			gamn= gam(ii)
			DO j= 1, nsa
				ah(j,i)= gamn*ah(j,i)
			END DO
			ah(i,i)= ah(i,i) + am(ii)
		END DO

		CALL FACTOR( nsa, ah(1:nsa,1:nsa), ifail)
        	IF (ifail .EQ. 1 ) CALL FAFERR(910)

!-----------------------------------Form the error in slip stress
       		DO  i= 1, nsa
			ii= isa(i)
			taufg= s1(ii) * SPOWER( ( gam(ii)/dt ), am(ii))
			IF ((tau(ii)*taufg) .GT. 0.) THEN
				errtau(i)= ALOG( tau(ii)/taufg )
			ELSE
				errtau(i)= 0.
			ENDIF
		END DO

!--------------------------------RMS value of slip stress error
		tauerr= 0.
		DO i=1,nsa
			tauerr=tauerr+errtau(i)*errtau(i)
		END DO
		tauerr=SQRT(tauerr)

!--------------------------------Correction to slip increment, cap
		CALL SOLVER( nsa, ah(1:nsa,1:nsa), errtau(1:nsa))
		errmax=0.
		DO i= 1, nsa
			IF ( ABS(errtau(i)) .GT. errmax ) errmax= ABS(errtau(i))
		END DO

		IF (errmax .GT. facmax) THEN
			scale= facmax/errmax
		ELSE
			scale= 1.
		ENDIF
!--------------------------------Update slip increment estimate
		DO i= 1, nsa
			ii=isa(i)
			gam(ii)= gam(ii)*EXP(scale*errtau(i))
		END DO

		l_iter= (iter .LT. maxit) .AND. (tauerr .GT. conv)
		
	END DO
	
    IF (tauerr .GT. 0.5)  CALL FAFERR(987)  !fail, too big an error

!--------------------Get and complete ah matrix in X visco-plastic form
	CALL PLASTC( dt, am)

    	DO i= 1, nsa
		ii= isa(i)
		DO j=1, nsa
			jj= isa(j)
			ah(j,i)= tau(jj)*ah(j,i)
		END DO

		IF ( ABS(gam(ii)) .GT. tiny ) THEN
			ah(i,i)=ah(i,i) + am(ii)*tau(ii)/gam(ii)
		ELSE
			ah(i,i)=ah(i,i) + am(ii)/tiny
		ENDIF

	END DO 

END SUBROUTINE SLPFTA

!***********************************************************************!
!																		!
!      SUBROUTINE SLPINI( dt, am)										!
!             															!
!-----------------------------------------------------------------------!
!      Routine for slip increments from tau for given resistances		!
!      with a cut-off for small relative slip increments to zero.		!
!-----------------------------------------------------------------------!
!   Arguments:	dt		REAL;				Time step					!
!               am		REAL;				Rate sensitivity index		!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE SLPINI(jgp, dt, am)

	USE FAS_COM			! FASOLT main variables declaration
	USE FAS_CPL			! Crystal plasticity variables

	IMPLICIT NONE

	REAL(4)		dt, twinrn

	INTEGER(4)	i, jgp
	REAL(4)		cofact, gammax, cutoff, SPOWER
	REAL(4), DIMENSION(1:nss(mat_code)):: am


    DATA cofact/1.E-6/					!	Relative slip cut-off to zero

    am(1:nss(mat_code))=ampar(1,1:nss(mat_code),mat_code)	!  Rate sensitivities

!--------------------------Get slip increments and maximum absolute value
    gammax= 0.
    DO i= 1, nss(mat_code)

		IF ( abs(ssang(i,mat_code)) .gt. 1.5 ) THEN   ! twin system, don't allow slip, just reorient!
			gam(i)= 0.
			! check if twin crit. fulfilled (RSS > CRSS) and (not yet reached num. limit) and (not twinning on other sys.):
			if ( ((tau(i) .gt. s0(i)) .and. (timestwinned(jgp) .ne. twinlim) .and. (twinsys(jgp) .eq. 0)) ) then
				CALL RANDOM_NUMBER(twinrn)
				if (twinrn .lt. twinprob) then
					tflag = 1
					twinsys(jgp) = i
				end if
			end if
		ELSE IF ( ssang(i,mat_code)*tau(i)/s0(i) .LT. 0. ) THEN		! if ssang=1, don't allow slip for negative tau/s0!
			gam(i)= 0.
		ELSE
			gam(i)=dt*SPOWER((tau(i)/s0(i)),(1./am(i)))
		ENDIF


		IF ( ABS( gam(i) ) .GT. gammax ) gammax=ABS(gam(i))
    END DO

!---------Cut off 'inactive' systems, form pointer and scale slips
    cutoff= cofact*gammax
    nsa= 0
    DO i= 1, nss(mat_code)
		IF ( ABS(gam(i)) .LT. cutoff ) THEN
			gam(i)= 0.
		ELSE
			nsa= nsa +1
			isa(nsa)= i
		ENDIF
	END DO
  
END SUBROUTINE SLPINI

!***********************************************************************!
!																		!
!       REAL FUNCTION SPOWER( a, b)										!
!-----------------------------------------------------------------------!
!    Function for signed power with traps								!
!-----------------------------------------------------------------------!
!   Arguments:       a     REAL ; argument								!
!                    b     REAL ; exponent								!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
REAL FUNCTION SPOWER( a, b)

	IMPLICIT NONE

	REAL(4)		a, b

	REAL(4)		aa, allm, bllm, sparg
	DATA allm, bllm/ 80.0, 1.E-35/

    aa=ABS(a)

    IF ( aa .LT. bllm ) THEN
		SPOWER=0.
    ELSE
		sparg= b* ALOG(aa)
		IF (sparg .LT. -allm) sparg= -allm
		IF (sparg .GT.  allm) sparg=  allm
		SPOWER=SIGN(EXP(sparg), a)
    ENDIF

END FUNCTION SPOWER

!***********************************************************************!
!																		!
!     SUBROUTINE FACTOR( n, w, ifail)									!
!																		!
!-----------------------------------------------------------------------!
!      -LU Factoring of square matrix, no frills!						!
!-----------------------------------------------------------------------!
!  Arguments:   n      INTEGER, Order of matrix					!
!               w      REAL ARRAY(1:n,1:n); Matrix, LU factored on exit	!
!               ifail  INTEGER, IF 0 OK ELSE indeterminate				!
!																		!
!***********************************************************************!
SUBROUTINE FACTOR( n, w, ifail)

	IMPLICIT NONE

	INTEGER(4)		n, ifail
	REAL(4)			w
	INTEGER(4)		i, j, k						
      
	DIMENSION w( 1:n, 1:n)

    	ifail=1
    	DO k= 1, n-1
		DO  i= k+1, n
			IF ( w(k,k) .EQ. 0. )      RETURN
			w(i,k)= w(i,k)/w(k,k)
			DO  j= k+1, n
				 w(i,j)=w(i,j) - w(i,k)*w(k,j)
			END DO
		END DO
	END DO

	IF ( w(n,n) .EQ. 0. )      RETURN

    ifail=0

END SUBROUTINE FACTOR

!***********************************************************************!
!																		!
!     SUBROUTINE SOLVER( n, w, x)										!
!																		!
!-----------------------------------------------------------------------!
!      -Linear equation solution given LU factored matrix				!
!-----------------------------------------------------------------------!
!  Arguments:   n      INTEGER, Order of matrix							!
!               w      REAL ARRAY(1:n,1:n); LU factored matrix			!
!               x      REAL ARRAY(n); LHS on entry, solution on exit	!
!																		!
!***********************************************************************!
SUBROUTINE SOLVER( n, w, x)

	IMPLICIT NONE

	INTEGER(4)	n
	REAL(4)		w, x

	INTEGER(4)	i, j

	DIMENSION	w(1:n,1:n),x(1:n)

	DO i= 2, n
		DO j= 1, i-1
			x(i)=x(i) - w(i,j)*x(j)
		END DO
	END DO

	x(n)=x(n)/w(n,n)

    	DO i= n-1, 1, -1
		DO j= i+1, n
			x(i)= x(i) - w(i,j)*x(j)
 		END DO
		x(i)= x(i)/w(i,i)
	END DO
      
END SUBROUTINE SOLVER

