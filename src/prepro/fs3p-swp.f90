!*******************************************************************************!
!																				!
!		SUBROUTINE SETSWTOL														!
!		 Sets sweeping tolerance to alpha* minimum intra-element distance.		!
!																				!
!*******************************************************************************!
SUBROUTINE SETSWTOL(tolerance)

	USE FAS_COM

	IMPLICIT NONE
	INTEGER(4)	i, j, k, l
	REAL(4)		tolerance, alpha, xlocal, lengths, del_x, lmin

	DIMENSION xlocal( 1:3, 1:20)

	DO	i= 1, nlmnt		!	loop on elements

		DO j= 1, nd(i)	!   get local coordinate list
			DO k=1, 3
				xlocal(k,j)= x( 3*ln(j,i) -3 + k )
			END DO
		END DO

		DO j= 1, nd(i)-1	! calculate distances, get minimum
			DO k= j+1, nd(i)

				lengths= 0.
				DO l= 1, 3
					del_x= xlocal(l,j)- xlocal(l,k)
					lengths= lengths + del_x*del_x
				END DO
				
				IF ( ((i+j+k) .EQ. 4 ) .OR. ( lengths .LT. lmin) )  lmin= lengths

			END DO
		END DO

	END DO

	alpha= 0.5
	tolerance= alpha *SQRT( lmin )	!reset alpha interactively at some stage


END SUBROUTINE

!*******************************************************************************!
!																				!
!		SUBROUTINE SWEEP														!
!		 Eliminates common nodes (within tolerance) from mesh					!
!																				!
!*******************************************************************************!
SUBROUTINE SWEEP( tolerance)

	USE FAS_COM

	IMPLICIT NONE
	LOGICAL(4)	lcommon
	INTEGER(4)	i_test, j, k, l
	REAL(4)		tolerance, lengths, del_x


	i_test= 1
	DO WHILE( i_test .LT. nnod )
		i_test= i_test + 1

		j=1				!	check with nodes up to less than one on trial
		lcommon=.FALSE.
		DO WHILE( (j .LT. i_test-1) .AND. (.NOT. lcommon) )
			j=j+1

			lengths= 0.
			DO k= -2, 0
				del_x= x(3*i_test + k) - x(3*j + k)
				lengths= lengths + del_x*del_x
			END DO

			lcommon=  SQRT(lengths) .LT. tolerance
			  
		END DO

			
		IF( lcommon) THEN  ! i_test is common, remove

			DO k= i_test, nnod-1	!	shuffle down node attributes
				DO l= -2, 0
					 x(3*k +l)=  x(3*(k+1) +l)
					dx(3*k +l)= dx(3*(k+1) +l)
					fc(3*k +l)= fc(3*(k+1) +l)
				END DO
				ib(k)= ib(k+1)
			END DO

			DO k= 1, nlmnt			!  re-allocate topology
				DO l= 1, nd(k)
					IF ( ln(l,k) .EQ. i_test ) ln(l,k)= j
					IF ( ln(l,k) .GT. i_test ) ln(l,k)= ln(l,k) -1
				END DO
			END DO

			nnod= nnod-1			!  reduce number of nodes
			i_test= i_test -1		!  set marker back one

		END IF

	END DO


END SUBROUTINE
