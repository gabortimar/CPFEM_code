!*******************************************************************************!
!																				!
!	SUBROUTINE SETBOUNDS														!
!		Simple brick boundary condition setting									!
!-------------------------------------------------------------------------------!
!  
!																				!
!*******************************************************************************!
SUBROUTINE SETBOUNDS(tolerance, xfree)

	USE FAS_COM

	IMPLICIT NONE
	
	INTEGER(4)	i, j, k, l_side, xfree
	REAL(4)		tolerance, xtreme
	DIMENSION	xtreme( 1:3, 1:2)

!-------------------------------------------------------------brick extreme coordinates
	xtreme(1:3,1)=x(1:3)
	xtreme(1:3,2)=x(1:3)

	DO i= 2, nnod
		DO j= 1, 3
			k= 3*i -3 +j
			IF ( x(k) .LT. xtreme(j,1) ) xtreme(j,1)=x(k)
			IF ( x(k) .GT. xtreme(j,2) ) xtreme(j,2)=x(k)
		END DO
	END DO

!-------------------------------------------------------------set BC codes on each node
	DO i= 1, nnod


!----------------------------------------------determine if at brick surface
		l_side= 0
		DO j= 1, 3
			k=3*i -3 +j		

			IF( ABS(x(k)-xtreme(j,1)) .LT. tolerance ) l_side= l_side + 2**( 2*j-2) 
			IF( ABS(x(k)-xtreme(j,2)) .LT. tolerance ) l_side= l_side + 2**( 2*j-1) 

		END DO

		IF ( l_side .GT. 0 ) THEN

!-----------------------------absolute displacements 'full' specifications as defaults
			SELECT CASE( l_side)   

				CASE( 1)	! x fixed (x min)
					ib(i)= 100
				CASE( 4)	! y fixed (y min)
					ib(i)= 200
				CASE(16)	! z fixed (z min)
					ib(i)= 300
				CASE( 5)	! x fix, y fix (x min, y min)
					ib(i)= 400
				CASE(17)	! x fix, z fix (x min, z min)
					ib(i)= 500
				CASE(20)	! y fix, z fix (y min, z min)
					ib(i)= 600
				CASE(21)	! all fixed (x min, y min, z min)
					ib(i)= 700

				CASE( 2)	! x spec. (x max)
					ib(i)= 800
				CASE( 8)	! y spec. (y max)
					ib(i)= 900
				CASE(32)	! z spec. (z max)
					ib(i)=1000
				CASE(10)	! x, y spec. (x max, y max)
					ib(i)=1100
				CASE(34)	! x, z spec. (x max, z max)
					ib(i)=1200
				CASE(40)	! y, z spec. (y max, z max)
					ib(i)=1300
				CASE(42)	! x, y, z spec. (x max, y max, z max)
					ib(i)=1400

				CASE( 9)	! x fix, y spec. (x min, y max)
					ib(i)=1500
				CASE(33)	! x fix, z spec. (x min, z max)
					ib(i)=1600
				CASE( 6)	! y fix, x spec. (y min, x max)
					ib(i)=1700
				CASE(36)	! y fix, z spec. (y min, z max)
					ib(i)=1800
				CASE(18)	! z fix, x spec. (z min, x max)
					ib(i)=1900
				CASE(24)	! z fix, y spec. (z min, y max)
					ib(i)=2000

				CASE(26)	! x, y spec., z fixed (x max, y max, z min)
					ib(i)=2100
				CASE(38)	! x, z spec., y fixed (x max, y min, z max)
					ib(i)=2200
				CASE(41)	! y, z spec., x fixed (x min, y max, z max)
					ib(i)=2300
				CASE(37)	! x, y fixed, z spec. (x min, y min, z max)
					ib(i)=2400
				CASE(25)	! x, z fixed, y spec. (x min, y max, z min)
					ib(i)=2500
				CASE(22)	! y, z fixed, x spec. (x max, y min, z min)
					ib(i)=2600

				CASE DEFAULT					!free boundary	(any true)
					ib(i)=1

			END SELECT

!--------------------------options for relative (strain), FC, 'channel die' or 'tensile' freedoms
			SELECT CASE( xfree )
				CASE(0)							! 'strain' condition
					 IF ( l_side .GT. 0) THEN
						ib(i)= 5000
					 END IF

				CASE(1)							! fully specified faces

				CASE(2)							! free x faces

					SELECT CASE (l_side)

						CASE( 1,2)
							ib(i)= 1	!  free ( x min or x max)

						CASE( 5,6)
							ib(i)= 200	!  y fixed ( x min, y min)
						CASE( 9,10)
							ib(i)= 900	!  y spec. ( x min, y max)

						CASE( 17,18)
							ib(i)= 300	!  z fixed ( x min, z min)
						CASE( 33,34)
							ib(i)=1000	!  z spec. ( x min, z max)

						CASE( 22)
							ib(i)= 600	! y fixed, z fixed (x max, y min, z min)
						CASE( 25,26)
							ib(i)=2000	! y spec., z fixed (x max, y max, z min)

						CASE( 37,38)
							ib(i)=1800	! y fixed, z spec. ( x min, y min or max, z max)
						CASE( 41,42)
							ib(i)=1300	! y spec., z spec. ( x min, y min or max, z max)

					END SELECT

				CASE(3)						! only x faces

					SELECT CASE (l_side)
						CASE(4,8,16,20,24,32,36,40)
							ib(i)=1		!	free ( not xmin or xmax)

						CASE(1,5,9,17,25,33,37,41)
							ib(i)=100	!   x fixed ( xmin)

						CASE(6,10,18,22,26,34,38,42)
							ib(i)=800	!	x specified ( xmax)


					END SELECT

				CASE DEFAULT

			END SELECT


		END IF

	END DO


END SUBROUTINE

!*******************************************************************************!
!																				!
!	SUBROUTINE SETKINC															!
!		Sets the  kinematic conditions											!
!																				!
!*******************************************************************************!
SUBROUTINE SETKINC( xl, xfree)

	USE FAS_COM

	IMPLICIT NONE
	INTEGER(4)	xfree, i, j
	REAL(4)		xl, eps_inc_X, eps_inc_Y, eps_inc_Z, q, xcurr, ycurr, zcurr, dxs, dys, dzs

	DIMENSION xl(1:3)

	SELECT CASE (xfree)
		CASE(0, 1)
			WRITE(*,'(/,5X,''Enter size of increments (e11, e22 and e33): '',\)')
			READ(*,*) eps_inc_X, eps_inc_Y, eps_inc_Z

		CASE(2)
			WRITE(*,'(/,5X,''Enter size of increments (e22 and e33): '',\)')
			READ(*,*) eps_inc_Y, eps_inc_Z
			eps_inc_X = 0.
			
		CASE(3) 
			WRITE(*,'(/,5X,''Enter size of increment (e11): '',\)')
			READ(*,*) eps_inc_X
			eps_inc_Y = 0.
			eps_inc_Z = 0.

	END SELECT

	WRITE(*,'(/,5X,''Enter number of increments (max. 2000): '',\)')
	READ(*,*) ninc1

!-------------------------------------------------boundary conditions
	nxs=50
	nfs=0

	x_spec=0.

	xcurr= xl(1)
	ycurr= xl(2)
	zcurr= xl(3)

	do i=1,ninc1
		tiem(i)= SQRT(eps_inc_X**2 + eps_inc_Y**2 + eps_inc_Z**2)	!	this gives mean rate approx. 1

		x_spec(1:6, 1,i)=(/1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
		x_spec(1:6, 2,i)=(/0.0, 0.0, 1.0, 0.0, 0.0, 0.0 /)
		x_spec(1:6, 3,i)=(/0.0, 0.0, 0.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6, 4,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 0.0 /)
		x_spec(1:6, 5,i)=(/1.0, 0.0, 0.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6, 6,i)=(/0.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6, 7,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)

		x_spec(1:6, 8,i)=(/1.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
		x_spec(1:6, 9,i)=(/0.0, 0.0, 1.0, 0.0, 0.0, 0.0 /)
		x_spec(1:6,10,i)=(/0.0, 0.0, 0.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6,11,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 0.0 /)
		x_spec(1:6,12,i)=(/1.0, 0.0, 0.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6,13,i)=(/0.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6,14,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)

		x_spec(1:6,15,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 0.0 /)
		x_spec(1:6,16,i)=(/1.0, 0.0, 0.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6,17,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 0.0 /)
		x_spec(1:6,18,i)=(/0.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6,19,i)=(/1.0, 0.0, 0.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6,20,i)=(/0.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)

		x_spec(1:6,21,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6,22,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6,23,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6,24,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6,25,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)
		x_spec(1:6,26,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)

		dxs = eps_inc_X*xcurr
		dys = eps_inc_Y*ycurr
		dzs = eps_inc_Z*zcurr

		xcurr= xcurr + dxs
		ycurr= ycurr + dys
		zcurr= zcurr + dzs

		DO	j= 1, 26
			x_spec(7:15,j,i)= 0.
		END DO

		x_spec( 7,  8, i)= dxs
		x_spec( 7, 11, i)= dxs
		x_spec( 7, 12, i)= dxs
		x_spec( 7, 14, i)= dxs
		x_spec( 7, 17, i)= dxs
		x_spec( 7, 19, i)= dxs
		x_spec( 7, 21, i)= dxs
		x_spec( 7, 22, i)= dxs
		x_spec( 7, 26, i)= dxs

		x_spec(11,  9, i)= dys
		x_spec(11, 11, i)= dys
		x_spec(11, 13, i)= dys
		x_spec(11, 14, i)= dys
		x_spec(11, 15, i)= dys
		x_spec(11, 20, i)= dys
		x_spec(11, 21, i)= dys
		x_spec(11, 23, i)= dys
		x_spec(11, 25, i)= dys


		x_spec(15, 10, i)= dzs
		x_spec(15, 12, i)= dzs
		x_spec(15, 13, i)= dzs
		x_spec(15, 14, i)= dzs
		x_spec(15, 16, i)= dzs
		x_spec(15, 18, i)= dzs
		x_spec(15, 22, i)= dzs
		x_spec(15, 23, i)= dzs
		x_spec(15, 24, i)= dzs

!-----------------------------------------full constraints specified straining
		x_spec(1:6,50,i)=(/1.0, 0.0, 1.0, 0.0, 0.0, 1.0 /)
		x_spec( 7, 50, i)= eps_inc_X
		x_spec( 9, 50, i)= 0.
		x_spec(11, 50, i)= eps_inc_Y
		x_spec(15, 50, i)= eps_inc_Z

	END DO

END SUBROUTINE

!*******************************************************************************!
!																				!
!	SUBROUTINE SETTIES															!
!		Includes tie elements for periodic boundary conditions					!
!																				!
!*******************************************************************************!
SUBROUTINE SETTIES( tolerance)

	USE FAS_COM

	IMPLICIT NONE
	
	CHARACTER(LEN=1) cdum
	LOGICAL(4)	tie
	INTEGER(4)	i, j, k, i1,i2,j1,j2,j3,k1,k2
	REAL(4)		tolerance, xtreme
	DIMENSION	tie(1:3), xtreme( 1:3, 1:2)

!-------------------------------------------------get conditions
	WRITE(*,'(//,10X,'' Set tied node conditions '')')

	WRITE(*,'(12X,'' Tie x faces (y/n) '',\)')
	READ(*,'(A)') cdum
	tie(1)=  cdum.EQ.'y' .OR. cdum.EQ.'Y'

	WRITE(*,'(12X,'' Tie y faces (y/n) '',\)')
	READ(*,'(A)') cdum
	tie(2)=  cdum.EQ.'y' .OR. cdum.EQ.'Y'

	WRITE(*,'(12X,'' Tie z faces (y/n) '',\)')
	READ(*,'(A)') cdum
	tie(3)=  cdum.EQ.'y' .OR. cdum.EQ.'Y'

!---------------------------brick extreme coordinates

	xtreme(1:3,1)=x(1:3)
	xtreme(1:3,2)=x(1:3)

	DO i= 2, nnod
		DO j= 1, 3
			k= 3*i -3 +j
			IF ( x(k) .LT. xtreme(j,1) ) xtreme(j,1)=x(k)
			IF ( x(k) .GT. xtreme(j,2) ) xtreme(j,2)=x(k)
		END DO
	END DO

!------------------------------------------set ties: first get lesser ends
	DO i1= 1, nnod

		DO j1= 1, 3
			k1=3*i1 -3 +j1	
			IF (tie(j1) .AND. (ABS(x(k1)-xtreme(j1,1)) .LT. tolerance) ) THEN
				j2= 1 + MOD(j1,3)
				j3= 1 + MOD(j2,3)

!-------------------------search for second ( greater) tie ends	
				DO i2= 1, nnod
					k2=3*i2 -3 + j1

					IF (	    (ABS(x(k2)-xtreme(j1,2)) .LT. tolerance) 		&
						& .AND.	(ABS(x(3*i1-3+j2)-x(3*i2-3+j2)) .LT. tolerance) &
						& .AND.	(ABS(x(3*i1-3+j3)-x(3*i2-3+j3)) .LT. tolerance) )  THEN
						
!------------------------include tie elements (code 1x, 2y, 3z, 4yz, 5zx, 6xy)

						SELECT CASE (ib(i1)/100)
							CASE(1,2,3)
								nlmnt= nlmnt+1; le(nlmnt)=0; lt(nlmnt)=j1+3
								nd(nlmnt)=2; np(nlmnt)=0
								ln(1,nlmnt)=i1; ln(2,nlmnt)=i2	

							CASE(6,18,20)
								nlmnt= nlmnt+1; le(nlmnt)=0; lt(nlmnt)= 1
								nd(nlmnt)=2; np(nlmnt)=0
								ln(1,nlmnt)=i1; ln(2,nlmnt)=i2
									
							CASE(5,16,19)
								nlmnt= nlmnt+1; le(nlmnt)=0; lt(nlmnt)= 2
								nd(nlmnt)=2; np(nlmnt)=0
								ln(1,nlmnt)=i1; ln(2,nlmnt)=i2	

							CASE(4,15,17)
								nlmnt= nlmnt+1; le(nlmnt)=0; lt(nlmnt)= 3
								nd(nlmnt)=2; np(nlmnt)=0
								ln(1,nlmnt)=i1; ln(2,nlmnt)=i2	

						END SELECT
								
					END IF
				END DO
			END IF
		END DO

	END DO

		


END SUBROUTINE
