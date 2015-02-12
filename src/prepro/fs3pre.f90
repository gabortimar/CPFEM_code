!*******************************************************************************!
!																				!
!							FS3PRE	Pete Bate									!
!			  Simple pre-processor for 3d fem development						!
!																				!
!*******************************************************************************!
PROGRAM FS3PRE

	USE FAS_COM

	IMPLICIT NONE
	CHARACTER	cdum
	INTEGER(4)	nein, ix, iy, iz, nlocs, arr_type, xfree, message
	REAL(4)		xl, x0, dxl, tolerance
	REAL clock_start, clock_finish
	DIMENSION	nein( 1:3), xl( 1:3), x0( 1:3), dxl( 1:3)

	INTEGER(4), ALLOCATABLE::  location(:,:)

!---------------------------------------get numbers and geometry
	WRITE(*,'(25X,''Basic FASOLT3 pre-processor'',//)')

	WRITE(*,'(20X,''Enter number of divisions in x: '',\)')
	READ(*,*) nein(1)
	WRITE(*,'(20X,''Enter number of divisions in y: '',\)')
	READ(*,*) nein(2)
	WRITE(*,'(20X,''Enter number of divisions in z: '',\)')
	READ(*,*) nein(3)
	WRITE(*,'(/)')

	WRITE(*,'(20X,''Enter x length: '',\)')
	READ(*,*) xl(1)
	WRITE(*,'(20X,''Enter y length: '',\)')
	READ(*,*) xl(2)
	WRITE(*,'(20X,''Enter z length: '',\)')
	READ(*,*) xl(3)



!-----------------------------element set type and estimate of number of elements
!								-location array gives position of 'grain'
!
	WRITE(*,'(/,15X,''Select                  8 node bricks    (1) or'')')
	WRITE(*,'(15X,''                       20 node bricks    (2) or'')')
	WRITE(*,'(15X,''rhomb. dodec. array of 10 node tetraheda (3) : '',\)')
	READ(*,*) arr_type
	CALL cpu_time (clock_start)

	SELECT CASE ( arr_type)
		CASE(1,2)
			nlocs= nein(1)*nein(2)*nein(3)
		CASE(3)
			nlocs= 12 + 17*(nein(1)*nein(2)*nein(3)/2)
		CASE DEFAULT
			STOP
	END SELECT

	ALLOCATE ( location( 1:3, 1:nlocs) )

!-----------------------------------------  initialise
	nnod=  0
	nlmnt= 0
	ngps=  0

!-----------------------------------------  add elements
	dxl(1:3)= xl(1:3)/ REAL(nein(1:3),4)

	DO iz= 1, nein(3)
		x0(3)= dxl(3)* REAL( iz-1, 4)

		DO iy= 1, nein(2)			! nein(i) is the number of divisions in the i-th direction
			x0(2)= dxl(2)* REAL( iy-1, 4)

			DO ix= 1, nein(1)
				x0(1)= dxl(1)* REAL( ix-1, 4)

				SELECT CASE ( arr_type)
					CASE(1)						! add a linear brick

						CALL ADDEL_8NQ( x0, dxl)

						location( 1, nlmnt)= ix			! indicial location of nlmnt-th element
						location( 2, nlmnt)= iy
						location( 3, nlmnt)= iz

					CASE(2)						! add a quadratic brick

						CALL ADDEL_20NQ( x0, dxl)

						location( 1, nlmnt)= ix
						location( 2, nlmnt)= iy
						location( 3, nlmnt)= iz

					CASE(3)						!add a r.d. group of quadratic tetraheda

						IF ( MOD((ix+iy+iz),2) .EQ. 1) THEN	! 'middle' type group

							CALL ADDEL_RDGM( x0, dxl)

							location( 1, nlmnt-4:nlmnt)= ix
							location( 2, nlmnt-4:nlmnt)= iy
							location( 3, nlmnt-4:nlmnt)= iz

						ELSE				! 'corners' type groups

							CALL ADDEL_RDGC( x0, dxl)

							location( 1, nlmnt-11:nlmnt-10)= ix-1
							location( 2, nlmnt-11:nlmnt-10)= iy
							location( 3, nlmnt-11:nlmnt-10)= iz

							location( 1, nlmnt -9:nlmnt -8)= ix
							location( 2, nlmnt -9:nlmnt -8)= iy-1
							location( 3, nlmnt -9:nlmnt -8)= iz

							location( 1, nlmnt -7:nlmnt -6)= ix
							location( 2, nlmnt -7:nlmnt -6)= iy
							location( 3, nlmnt -7:nlmnt -6)= iz-1

							location( 1, nlmnt -5:nlmnt -4)= ix+1
							location( 2, nlmnt -5:nlmnt -4)= iy
							location( 3, nlmnt -5:nlmnt -4)= iz

							location( 1, nlmnt -3:nlmnt -2)= ix
							location( 2, nlmnt -3:nlmnt -2)= iy+1
							location( 3, nlmnt -3:nlmnt -2)= iz

							location( 1, nlmnt -1:nlmnt   )= ix
							location( 2, nlmnt -1:nlmnt   )= iy
							location( 3, nlmnt -1:nlmnt   )= iz+1

						END IF

				END SELECT

				IF (nnod .GE. maxnod-120) THEN	!	sweep now!
					CALL SETSWTOL(tolerance)	!	node sweep tolerance
					CALL SWEEP(tolerance)		!	sweep common nodes
				END IF

			END DO
		END DO
	END DO

	CALL SETSWTOL(tolerance)			!	node sweep tolerance
	CALL SWEEP(tolerance)				!	sweep common nodes

!---------------------------------------report number of nodes and elements


	WRITE(*,'(/,20X,''Number of nodes:'',I10,'' and elements: '',I5)') nnod, nlmnt
	call cpu_time(clock_finish)
	write (*,'(A)') "--------------------------------------"
	write (*,'(A)') "  "
	write (*,*) "CPU Time = ",(clock_finish - clock_start), " seconds"

!-----------------------------------------------boundary conditions
	WRITE(*,'(/,10X, ''Boundary strain(0) FC brick(1), Free x faces(2) or Just x faces(3)?'',\)')
	READ(*,*) xfree
	CALL SETBOUNDS(tolerance, xfree)

	CALL SETSTATE( nein, location )			!	set element state

	CALL SETMATC							!	set material data

	CALL SETKINC( xl, xfree)				!	set displacement conditions
	! read twinlim and twinprob variables:
	READ(*,*) twinlim
	READ(*,*) twinprob

	CALL SETTIES(tolerance)					!   include 'tie' elements

	CALL WTDATA					!	write file

END

!*******************************************************************************!
!																				!
!	SUBROUTINE ADDEL_8NQ														!
!		Adds in a 8 node hexahedron, orthogonal edges							!
!																				!
!*******************************************************************************!
SUBROUTINE ADDEL_8NQ( x0, delx )

	USE FAS_COM

	IMPLICIT NONE
	INTEGER(4)	i, j
	REAL(4)		x0, delx, xyz
	DIMENSION	x0( 1:3), delx( 1:3), xyz( 1:3, 1:8)

!-------------------------------------------------data: local coordinates
	DATA ( ( xyz(j,i),j=1,3),i=1,8)/  0.0, 0.0, 0.0,   1.0, 0.0, 0.0,	&
								&     1.0, 1.0, 0.0,   0.0, 1.0, 0.0,   &
								&     0.0, 0.0, 1.0,   1.0, 0.0, 1.0,	&
								&	  1.0, 1.0, 1.0,   0.0, 1.0, 1.0    /

	nlmnt= nlmnt +1
	ngps=  ngps  +1

	IF( nnod+8.GT.maxnod .OR. nlmnt.GT.maxels .OR. ngps.GT.maxgps ) THEN
			PRINT*,'   Too big!'
			PRINT*, nnod, nlmnt, ngps
			PAUSE
			STOP
	END IF


	le( nlmnt)= 1			
	lt( nlmnt)= 12 
	nd( nlmnt)= 8
	np( nlmnt)= 1



	DO i= 1, 8					!	nodes, initialise dx,fc and ib
		nnod= nnod+1

		ln( i, nlmnt)= nnod
		ib( nnod)= 0

		DO j= 1, 3
			x( 3*nnod -3 +j ) = x0(j) + delx(j)*xyz(j,i)

			dx( 3*nnod -3 +j) = 0.
			fc( 3*nnod -3 +j) = 0.
		END DO

	END DO
	
	lg( 1, nlmnt) = nlmnt		!	gauss point

END SUBROUTINE

!*******************************************************************************!
!																				!
!	SUBROUTINE ADDEL_20NQ														!
!		Adds in a 20 node hexahedron, orthogonal edges							!
!																				!
!*******************************************************************************!
SUBROUTINE ADDEL_20NQ( x0, delx )

	USE FAS_COM

	IMPLICIT NONE
	INTEGER(4)	i, j
	REAL(4)		x0, delx, xyz
	DIMENSION	x0( 1:3), delx( 1:3), xyz( 1:3, 1:20)

!-------------------------------------------------data: local coordinates
	DATA ( ( xyz(j,i),j=1,3),i=1,20)/  0.0, 0.0, 0.0,  0.5, 0.0, 0.0,  1.0, 0.0, 0.0,	&
	&  1.0, 0.5, 0.0,  1.0, 1.0, 0.0,  0.5, 1.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.5, 0.0,	&
	&  0.0, 0.0, 0.5,  1.0, 0.0, 0.5,  1.0, 1.0, 0.5,  0.0, 1.0, 0.5,  0.0, 0.0, 1.0,   &
	&  0.5, 0.0, 1.0,  1.0, 0.0, 1.0,  1.0, 0.5, 1.0,  1.0, 1.0, 1.0,  0.5, 1.0, 1.0,	&
	&  0.0, 1.0, 1.0,  0.0, 0.5, 1.0 /

	nlmnt= nlmnt +1		! current number of added elements
	ngps=  ngps  +8		! current number of added GPs

	IF( nnod+20.GT.maxnod .OR. nlmnt.GT.maxels .OR. ngps.GT.maxgps ) THEN
			PRINT*,'   Too big!'
			PRINT*, nnod, nlmnt, ngps
			PAUSE
			STOP
	END IF


	le( nlmnt)= 1			
	lt( nlmnt)= 22 
	nd( nlmnt)= 20		! number of nodes for current element
	np( nlmnt)= 8		! number of GPs for current element



	DO i= 1, 20					!	nodes, initialise dx,fc and ib
		nnod= nnod+1

		ln( i, nlmnt)= nnod			!	index of i-th node of nlmnt-th element
		ib( nnod)= 0

		DO j= 1, 3
			x( 3*nnod -3 +j ) = x0(j) + delx(j)*xyz(j,i)	! set node position

			dx( 3*nnod -3 +j) = 0.
			fc( 3*nnod -3 +j) = 0.
		END DO

	END DO
	
	DO i= 1, 8					!	gauss points
		lg( i, nlmnt) = 8*nlmnt -8 +i		!	index of i-th GP of nlmnt-th element
	END DO	

END SUBROUTINE

!*******************************************************************************!
!																				!
!	SUBROUTINE ADDEL_RDGM														!
!  Adds in a middle group of a rhombic dodecahedral array of 10 node tetrahedra	!
!																				!
!*******************************************************************************!
SUBROUTINE ADDEL_RDGM( x0, delx )

	USE FAS_COM

	IMPLICIT NONE
	INTEGER(4)	i, j, ie
	REAL(4)		x0, delx, xt
	DIMENSION	x0( 1:3), delx( 1:3), xt( 1:3, 1:10, 1:5)

!---------------------------------------------data: local coordinates
!-----------------------------------------------------------------------------------!
!    Data for the type 1 discretisation. The elements locations are:				!
!      1 middle; 2 all low; { 3 x low, 4 y low, 5 z low (rest high)	}				!
!-----------------------------------------------------------------------------------!

	DATA ((xt(i,j,1),i=1,3),j=1,10)/	&
	&	0.00,1.00,0.00,0.00,0.50,0.50,0.00,0.00,1.00,0.50,0.00,0.50,1.00,0.00,0.00,	&
	&	0.50,0.50,0.00,0.50,1.00,0.50,0.50,0.50,1.00,1.00,0.50,0.50,1.00,1.00,1.00	/

	DATA ((xt(i,j,2),i=1,3),j=1,10)/	&
	&	0.00,1.00,0.00,0.50,0.50,0.00,1.00,0.00,0.00,0.50,0.00,0.50,0.00,0.00,1.00,	&
	&	0.00,0.50,0.50,0.00,0.50,0.00,0.50,0.00,0.00,0.00,0.00,0.50,0.00,0.00,0.00	/

	DATA ((xt(i,j,3),i=1,3),j=1,10)/	&
	&	0.00,1.00,0.00,0.00,0.50,0.50,0.00,0.00,1.00,0.50,0.50,1.00,1.00,1.00,1.00,	&
	&	0.50,1.00,0.50,0.00,1.00,0.50,0.00,0.50,1.00,0.50,1.00,1.00,0.00,1.00,1.00	/

	DATA ((xt(i,j,4),i=1,3),j=1,10)/	&
	&	1.00,0.00,0.00,1.00,0.50,0.50,1.00,1.00,1.00,0.50,0.50,1.00,0.00,0.00,1.00,	&
	&	0.50,0.00,0.50,1.00,0.00,0.50,1.00,0.50,1.00,0.50,0.00,1.00,1.00,0.00,1.00	/

	DATA ((xt(i,j,5),i=1,3),j=1,10)/	&
	&	0.00,1.00,0.00,0.50,1.00,0.50,1.00,1.00,1.00,1.00,0.50,0.50,1.00,0.00,0.00,	&
	&	0.50,0.50,0.00,0.50,1.00,0.00,1.00,1.00,0.50,1.00,0.50,0.00,1.00,1.00,0.00	/

!-----------------------------------add five tetrahedra
	DO ie= 1, 5
		nlmnt= nlmnt +1
		ngps=  ngps  +4

		IF( nnod+10.GT.maxnod .OR. nlmnt.GT.maxels .OR. ngps.GT.maxgps ) THEN
			PRINT*,'   Too big!'
			PRINT*, nnod, nlmnt, ngps
			PAUSE
			STOP
		END IF


		le( nlmnt)= 1			
		lt( nlmnt)= 21
		nd( nlmnt)= 10
		np( nlmnt)= 4

		DO i= 1, 10					!	nodes, initialise dx,fc and ib
			nnod= nnod+1

			ln( i, nlmnt)= nnod
			ib( nnod)= 0

			DO j= 1, 3
				x( 3*nnod -3 +j ) = x0(j) + delx(j)*xt(j,i,ie)

				dx( 3*nnod -3 +j) = 0.
				fc( 3*nnod -3 +j) = 0.
			END DO

		END DO
	
		DO i= 1, 4					!	gauss points
			lg( i, nlmnt) = 4*nlmnt -4 +i
		END DO
			
	END DO

END SUBROUTINE

!*******************************************************************************!
!																				!
!	SUBROUTINE ADDEL_RDGC														!
!  Adds in a corner group of a rhombic dodecahedral array of 10 node tetrahedra	!
!																				!
!*******************************************************************************!
SUBROUTINE ADDEL_RDGC( x0, delx )

	USE FAS_COM

	IMPLICIT NONE
	INTEGER(4)	i, j, ie
	REAL(4)		x0, delx, xt
	DIMENSION	x0( 1:3), delx( 1:3), xt( 1:3, 1:10, 1:12)

!-------------------------------------------------data: local coordinates
!-----------------------------------------------------------------------------------!
!    Data for the type 2 discretisation. The faces which the elements adjoin are:	!
!      1,2 low x; 3,4 low y; 5,6 low z;  7,8 high x;  9,10 high y;  11,12 high z	!
!-----------------------------------------------------------------------------------!

DATA ((xt(i,j, 1),i=1,3),j=1,10)/	&
	&	0.00,0.00,0.00,0.00,0.50,0.00,0.00,1.00,0.00,0.00,1.00,0.50,0.00,1.00,1.00,	&
	&	0.00,0.50,0.50,0.25,0.25,0.25,0.25,0.75,0.25,0.25,0.75,0.75,0.50,0.50,0.50	/

DATA ((xt(i,j, 2),i=1,3),j=1,10)/	&
	&	0.00,0.00,0.00,0.00,0.50,0.50,0.00,1.00,1.00,0.00,0.50,1.00,0.00,0.00,1.00,	&
	&	0.00,0.00,0.50,0.25,0.25,0.25,0.25,0.75,0.75,0.25,0.25,0.75,0.50,0.50,0.50	/

DATA ((xt(i,j, 3),i=1,3),j=1,10)/	&
	&	0.00,0.00,0.00,0.00,0.00,0.50,0.00,0.00,1.00,0.50,0.00,1.00,1.00,0.00,1.00,	&
	&	0.50,0.00,0.50,0.25,0.25,0.25,0.25,0.25,0.75,0.75,0.25,0.75,0.50,0.50,0.50	/

DATA ((xt(i,j, 4),i=1,3),j=1,10)/	&
	&	0.00,0.00,0.00,0.50,0.00,0.50,1.00,0.00,1.00,1.00,0.00,0.50,1.00,0.00,0.00,	&
	&0.50,0.00,0.00,0.25,0.25,0.25,0.75,0.25,0.75,0.75,0.25,0.25,0.50,0.50,0.50	/

DATA ((xt(i,j, 5),i=1,3),j=1,10)/	&
	&	0.00,0.00,0.00,0.50,0.00,0.00,1.00,0.00,0.00,1.00,0.50,0.00,1.00,1.00,0.00,	&
	&	0.50,0.50,0.00,0.25,0.25,0.25,0.75,0.25,0.25,0.75,0.75,0.25,0.50,0.50,0.50	/

DATA ((xt(i,j, 6),i=1,3),j=1,10)/	&
	&	0.00,0.00,0.00,0.50,0.50,0.00,1.00,1.00,0.00,0.50,1.00,0.00,0.00,1.00,0.00,	&
	&	0.00,0.50,0.00,0.25,0.25,0.25,0.75,0.75,0.25,0.25,0.75,0.25,0.50,0.50,0.50	/

DATA ((xt(i,j, 7),i=1,3),j=1,10)/	&
	&	1.00,1.00,0.00,1.00,0.50,0.00,1.00,0.00,0.00,1.00,0.00,0.50,1.00,0.00,1.00,	&
	&	1.00,0.50,0.50,0.75,0.75,0.25,0.75,0.25,0.25,0.75,0.25,0.75,0.50,0.50,0.50	/

DATA ((xt(i,j, 8),i=1,3),j=1,10)/	&
	&	1.00,1.00,0.00,1.00,0.50,0.50,1.00,0.00,1.00,1.00,0.50,1.00,1.00,1.00,1.00,	&
	&	1.00,1.00,0.50,0.75,0.75,0.25,0.75,0.25,0.75,0.75,0.75,0.75,0.50,0.50,0.50	/

DATA ((xt(i,j, 9),i=1,3),j=1,10)/	&
	&	1.00,1.00,0.00,1.00,1.00,0.50,1.00,1.00,1.00,0.50,1.00,1.00,0.00,1.00,1.00,	&
	&	0.50,1.00,0.50,0.75,0.75,0.25,0.75,0.75,0.75,0.25,0.75,0.75,0.50,0.50,0.50	/

DATA ((xt(i,j,10),i=1,3),j=1,10)/	&
	&	1.00,1.00,0.00,0.50,1.00,0.50,0.00,1.00,1.00,0.00,1.00,0.50,0.00,1.00,0.00,	&
	&	0.50,1.00,0.00,0.75,0.75,0.25,0.25,0.75,0.75,0.25,0.75,0.25,0.50,0.50,0.50	/

DATA ((xt(i,j,11),i=1,3),j=1,10)/	&
	&	1.00,0.00,1.00,0.50,0.00,1.00,0.00,0.00,1.00,0.00,0.50,1.00,0.00,1.00,1.00,	&
	&	0.50,0.50,1.00,0.75,0.25,0.75,0.25,0.25,0.75,0.25,0.75,0.75,0.50,0.50,0.50	/

DATA ((xt(i,j,12),i=1,3),j=1,10)/	&
	&	1.00,0.00,1.00,0.50,0.50,1.00,0.00,1.00,1.00,0.50,1.00,1.00,1.00,1.00,1.00,	&
	&	1.00,0.50,1.00,0.75,0.25,0.75,0.25,0.75,0.75,0.75,0.75,0.75,0.50,0.50,0.50	/

!-----------------------------------add twelve tetrahedra
	DO ie= 1, 12		

		nlmnt= nlmnt +1
		ngps=  ngps  +4

		IF( nnod+10.GT.maxnod .OR. nlmnt.GT.maxels .OR. ngps.GT.maxgps ) THEN
			PRINT*,'   Too big!'
			PRINT*, nnod, nlmnt, ngps
			PAUSE
			STOP
		END IF


		le( nlmnt)= 1
		lt( nlmnt)= 21
		nd( nlmnt)= 10
		np( nlmnt)= 4

		DO i= 1, 10					!	nodes, initialise dx,fc and ib
			nnod= nnod+1

			ln( i, nlmnt)= nnod
			ib( nnod)= 0

			DO j= 1, 3
				x( 3*nnod -3 +j ) = x0(j) + delx(j)*xt(j,i,ie)

				dx( 3*nnod -3 +j) = 0.
				fc( 3*nnod -3 +j) = 0.
			END DO

		END DO
	
		DO i= 1, 4					!	gauss points
			lg( i, nlmnt) = 4*nlmnt -4 +i
		END DO
			
	END DO
END SUBROUTINE

