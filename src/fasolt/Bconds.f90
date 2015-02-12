!***********************************************************************!
!																		!
!       SUBROUTINE SUVALS												!
!-----------------------------------------------------------------------!
!       Evaluation and assembly of boundary conditions and surface		!
!		 related quantities for mechanical problem.						!
!-----------------------------------------------------------------------!
!    BC codes:  ib code formulated as follows:							!
!				decimal ccbbaa											!
!				a: 0 nul, 1 surface node specifier						!
!				b: 0 nul, 1....99  Kinematic constraints				!
!				c: 0 nul, 1....99  Specified Loads						!
!																		!
!    Notes:																!
!			The kinematic stiffness is added: this allows freedoms      !
!		            (9, i_kin, inc): 1,6 'stiffness', 7..9 spec. dx		!
!																		!
!***********************************************************************!
SUBROUTINE SUVALS

	USE FAS_COM			! FASOLT main variables declaration

	IMPLICIT NONE

	INTEGER(4)	i, j, k, l, ndof, idof, inod, ibc,	&
		&		i_kin, i_spf, row, loctn, point
	REAL(4)		xsloc, xspec, d_spec, sum

	DIMENSION	point( 1:3, 1:3), xsloc(1:15)
	DIMENSION	xspec( 1:3), d_spec(1:3,1:3)

!--------------------------------------Pointer to UT 3x3 array
	DATA ((point(i,j), i=1,3), j=1,3)/1,2,4, 2,3,5, 4,5,6/

!------------------------------Put current inbalances to r.h.s. vectors
	ndof= 3*nnod
    ddx(1:ndof)=fc(1:ndof)


	DO  inod= 1, nnod
		ibc= ib(inod)
		IF ( ibc .GE. 100 ) THEN		! There is a boundary condition
			idof= 3*inod-3

			i_kin= MOD( ibc/100, 100)
			IF ( i_kin .GT. 0 )	THEN		!	Kinematic constraint
				
				xsloc(1:15)= x_spec(1:15, i_kin, inc)

				IF ( i_kin.LT. 50) THEN		! absolute constraint
					xspec(1)= xsloc( 7)
					xspec(2)= xsloc(11)
					xspec(3)= xsloc(15)

				ELSE						! relative constraint	
					l= 6
					DO j= 1,3
						DO k=1,3
						  l=l+1
						  d_spec(j,k)= xsloc(l)
						END DO
					END DO

					DO j= 1, 3
						xspec(j)=0.
						DO k= 1, 3
							xspec(j)= xspec(j) + d_spec(j,k)*( x(idof+k) + dx(idof+k) )
						END DO
					END DO

				END IF

!--------------------------------------augment stiffness, apply displacement
				DO i= 1, 3

					row= idof + i
					DO j= 1, 3
						IF ( j .NE. i ) THEN
							k=2
							DO WHILE( idof+j .NE. is(k, row) )
								k=k+1
							END DO
						ELSE
							k=1
						END IF
						ss(k, row)= ss(k, row) +(efix*xsloc(point(i,j)))
					END DO

					sum=0.
					DO j= 1, 3
						sum=sum + efix* xsloc(point(i,j)) * ( dx(idof+j) -xspec(j) )
					END DO

					ddx( row)= ddx( row) + sum

				END DO

			END IF

			i_spf= MOD( ibc/10000, 100)
			IF ( i_spf .GT. 0 )	THEN		!	Specified load

				DO i= 1, 3
					ddx( idof+i)= ddx( idof+i) + f_spec( i, i_spf, inc)
				END DO

			END IF

		END IF
	END DO


END SUBROUTINE SUVALS