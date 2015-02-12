!***********************************************************************!
!																		!
!          SUBROUTINE NRITER											!
!-----------------------------------------------------------------------!
!			Damped Newton-Raphson iteration for solution				!
!			of deformation over a single increment.						!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE NRITER

	USE FAS_COM			! FASOLT main variables declaration
	USE FAS_CPL			! Crystal plasticity variables

	IMPLICIT NONE

	LOGICAL(4)		l_conv
	INTEGER(4)		nit, ndof, i, igp, iop
	REAL(4)			sfc, erddx, errdx, alp2, xc, sdx
     
!---------------------------------------Common for sum of iteration error
	REAL(4)		sumcpiter
	COMMON/ CPITER / sumcpiter

!--------------------write header to monitor report
	CALL ITHEAD(inc)
!--------------------------Initiate repeat sequence
	ndof=3*nnod
    nit= 0
	l_conv= .TRUE.

	DO WHILE(l_conv)
		nit= nit + 1

		sumcpiter= 0.	! initialise CP iteration error sum

!-------------------Initialise mechanical arrays
        fc(1:ndof)= 0.
        ss( 1:1000, 1:ndof)= 0.

!------------------Form equations for the corrections
        CALL ELVALS
        CALL SUVALS

!------------------Get internal imbalances of F
		sfc=0.
		DO i= 1, nnod
			IF ( ib(i).LT.10 ) 	sfc=sfc + ddx(3*i-2)*ddx(3*i-2) +	&
							&	  ddx(3*i-1)*ddx(3*i-1) + ddx(3*i)*ddx(3*i)
		END DO
		sfc= SQRT(sfc)
!-----------------------------------------------Solve
		CALL FACSLS(ndof)

!-------Calculate displacement error and damping parameter
		errdx=SQRT( DOT_PRODUCT( dx(1:ndof), dx(1:ndof)) )
		erddx=SQRT( DOT_PRODUCT(ddx(1:ndof),ddx(1:ndof)) )
		IF( errdx .EQ. 0. )  THEN
			sdx= toler
		ELSE
			sdx= erddx/errdx
		ENDIF
		alp2=1./(1.+alpha*sdx)

		sumcpiter= sumcpiter/ REAL( ngps, 4)

		CALL ITEREP(nit,sdx,sfc,sumcpiter,alp2)	!iteration progress report

!------------------------------Correct incremental values
		DO i= 1, ndof
			xc= alp2*ddx(i)
			x(i)=  x(i)  - xc
			dx(i)= dx(i) - xc
		END DO

       l_conv= (sdx.GE.toler).AND.(nit.LT.maxit)
	END DO

	! update cumslip variable here:
	do igp=1, ngps
	  do iop=1,nss(mat_code)
		cumslip(iop, igp) = cumslip(iop, igp) + (srate(iop, igp) * tiem(inc))	! slip added up
	  end do
	end do

	CALL POWER		! report power

END SUBROUTINE NRITER

!***********************************************************************!
!																		!
!          SUBROUTINE SETPTA											!
!-----------------------------------------------------------------------!
!	Generation of pointer arrays for compact storage of stiffness		!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE SETPTA

	USE FAS_COM			! FASOLT main variables declaration

	IMPLICIT NONE

	LOGICAL(4)	new
	INTEGER(4)	ilmnt, i, j, k, js, ks, locpoint, row, column, noffdiag
	DIMENSION	locpoint(1: maxofd)

	is(1:maxofd,1:3*nnod) = 0		! initialise pointer

	DO ilmnt= 1, nlmnt			! loop on elements

		DO i= 1, nd(ilmnt)		! loop on its nodes, form DOF array
			DO j= -2, 0 
				locpoint( 3*i + j) = 3*ln(i,ilmnt) + j
			END DO
		END DO

		DO	i= 1, 3*nd(ilmnt)				! loop on pointer row dof
			row= locpoint(i)

			DO j= 1, 3*nd(ilmnt)
				column=locpoint(j)

			    IF( row .NE. column )  THEN		! loop, off-diagonal column dofs

					new=.TRUE.					! check if already in list
					DO k= 2, 1 + is( 1, row  )
						IF( is( k, row) .EQ. column) new=.FALSE.
					END DO

					IF ( new ) THEN				! not in list, add
						noffdiag= is(1,row) + 1

							IF (noffdiag .GT. maxofd)  CALL FAFERR(955)

						is( 1, row)= noffdiag
						is( 1+noffdiag, row)= column
					END IF

				END IF

			END DO


		END DO

	END DO

END SUBROUTINE SETPTA

!***********************************************************************!
!																		!
!         SUBROUTINE FACSLS(M)											!
!-----------------------------------------------------------------------!
!				Solution of linear equations using SOR					!
!-----------------------------------------------------------------------!
!         ARGUMENTS           m     INTEGER, SIZE OF R.H.S VECTOR		!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE FACSLS(m)

	USE FAS_COM			! FASOLT main variables declaration

	IMPLICIT NONE

	LOGICAL(4)	conv
	INTEGER(4)	m, iter, maxitl, i, j 
	REAL(4)		tolerl, alphal, y, dxl, error
	REAL(4)		sum
	DIMENSION	y(1:maxdof), dxl(1:maxdof)

	DATA	alphal, maxitl/ 1.2, 1000/

	tolerl= toler*toler

!---------allocate lhs and rhs vectors, set convergence flag
	y(1:m)= ddx(1:m)
	ddx(1:m)= 0.
	conv= .FALSE.

!---------------------------------------------test diagonals
	DO i= 1, m
		IF ( ABS( ss(1,i)) .EQ. 0.) CALL FAFERR(903)
	END DO
	
!----------------------------------------------------iterate
	iter=0
	DO WHILE( (.NOT. conv) .AND. (iter.LT.maxitl))
		iter= iter + 1
		dxl(1:m)= ddx(1:m)				! store old x

		DO i= 1, m						! sweep on rows

!------------go sideways to get off-diagonal a.x terms
			sum=0.D0
			DO j=2, 1 + is(1,i)
				sum= sum + ss(j,i)*ddx(is(j,i))
			END DO

!------------------get change term, add weighted value
			sum=( y(i) - sum )/ss(1,i)
			ddx(i)= ddx(i)+ alphal*( sum - ddx(i))

		END DO

!--------------------get change in x, test convergence
		dxl= ddx- dxl
		error= DOT_PRODUCT(dxl,dxl)/DOT_PRODUCT(ddx,ddx)
		conv= ( error .LT. tolerl )	
	END DO
	
	IF ( error .GT. 10.*tolerl)  CALL FAFERR(904)	! convergence fail	

END SUBROUTINE FACSLS

