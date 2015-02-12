!***********************************************************************************!
!																					!
!	REAL(4) FUNCTION MISORI															!
!-----------------------------------------------------------------------------------!
!			Misorientation angle		 						 					!
!																					!
!***********************************************************************************!
REAL(4) FUNCTION MISORI( isym, ea_ref, ea_val )

	IMPLICIT NONE
	INTEGER(4)	isym
	REAL(4)		ea_ref, ea_val, a, b, trmax, QTNTCE
	DIMENSION	ea_ref(1:3), ea_val(1:3), a(0:3), b(0:3)
	
	CALL QTNEUL( ea_ref(1), ea_ref(2), ea_ref(3), a)
	CALL QTNEUL( ea_val(1), ea_val(2), ea_val(3), b)

	CALL SYM_NEAR( isym, a, b, trmax)

    MISORI= 2.*ACOS( trmax)

END FUNCTION

!*******************************************************************************!
!																				!
!	SUBROUTINE QTNEUL(ph1,phi,ph2,x)											!
!-------------------------------------------------------------------------------!
!   Returns the Euler-symmetrics (0:3) given the Euler angles					!
!																				!
!*******************************************************************************!
SUBROUTINE QTNEUL(ph1,phi,ph2,x)

	IMPLICIT NONE
	REAL(4)		ph1, phi, ph2, x, c_half, s_half
	DIMENSION	x(0:3)

	c_half= COS( phi/2.)
	s_half= SIN( phi/2.)
	
	x(0)=  c_half * COS( (ph1+ph2)/2. )
	x(1)= -s_half * COS( (ph2-ph1)/2. )
	x(2)=  s_half * SIN( (ph2-ph1)/2. )
	x(3)= -c_half * SIN( (ph1+ph2)/2. )

	IF ( x(0) .LT. 0. ) x= -x		! ensure in positive hemisphere
	

END SUBROUTINE QTNEUL

!*******************************************************************************!
!																				!
!	SUBROUTINE EULQTN(a,ph1,phi,ph2)											!
!-------------------------------------------------------------------------------!
!   Returns the Euler angles given the Euler symmetrics, a(0:3)					!
!																				!
!*******************************************************************************!
SUBROUTINE EULQTN(a,ph1,phi,ph2)

	IMPLICIT NONE	
	REAL(4)		a, ph1, phi, ph2
	DIMENSION	a(0:3)

	REAL(4)	cphi

	cphi= 1. -2.*(a(1)*a(1)+a(2)*a(2))
	IF ( cphi .GT. -1.) THEN 
		phi= ACOS( cphi)
	ELSE
		phi= 3.14159265359
	END IF

	IF ( ABS(cphi) .GE. 1 ) THEN
		ph1= 0.5*ATAN2(  2.*(a(1)*a(2)-a(0)*a(3)) ,1.-2.*(a(2)*a(2)+a(3)*a(3)))
		ph2= ph1
	ELSE
		ph1= ATAN2( (a(1)*a(3)-a(0)*a(2)), -(a(2)*a(3)+a(0)*a(1)) )
		ph2= ATAN2( (a(1)*a(3)+a(0)*a(2)),  (a(2)*a(3)-a(0)*a(1)) )
	END IF
	
END SUBROUTINE EULQTN

!*******************************************************************************!
!																				!
!	REAL FUNCTION QTNTCE(a,b)													!
!-------------------------------------------------------------------------------!
!		Returns the scalar product two Euler-symmetrics a(0:3), b(0:3). 		!
!			-note the relationship with the first component of a':b.			!
!			( see QTNPRD and QTNTRP code )										!
!																				!
!*******************************************************************************!
REAL FUNCTION QTNTCE(a,b)

	IMPLICIT NONE
	INTEGER(4)	i
	REAL(4)		a, b, sum
	DIMENSION	a(0:3), b(0:3)

	sum= 0.
	DO i= 0, 3
		sum= sum + a(i)*b(i)
	END DO

	QTNTCE= sum

END FUNCTION QTNTCE

!*******************************************************************************!
!																				!
!	SUBROUTINE QTNNRM(a, norm)													!
!-------------------------------------------------------------------------------!
!   Returns the normalised Euler symmetrics (0:3) and the normalisation factor.	!
!   Note that the normalisation can be used to give 'inverse variance'			!
!																				!
!*******************************************************************************!
SUBROUTINE QTNNRM(a, norm)

	IMPLICIT NONE
	REAL(4)		a, norm, QTNTCE
	DIMENSION	a(0:3)

	norm= 1./ SQRT( QTNTCE(a,a) )

	a= norm * a

END SUBROUTINE QTNNRM

!*******************************************************************************!
!																				!
!	SUBROUTINE QTNPRD(a,b,c)													!
!-------------------------------------------------------------------------------!
!		Forms the product c = a:b, in Euler-symmetrics (0:3)					!
!																				!
!*******************************************************************************!
SUBROUTINE QTNPRD(a,b,c)

	IMPLICIT NONE
	REAL(4)		a, b, c
	DIMENSION	a(0:3), b(0:3), c(0:3)

	c(0)= a(0)*b(0) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
	c(1)= a(0)*b(1) + a(1)*b(0) + a(2)*b(3) - a(3)*b(2)
	c(2)= a(0)*b(2) - a(1)*b(3) + a(2)*b(0) + a(3)*b(1)
	c(3)= a(0)*b(3) + a(1)*b(2) - a(2)*b(1) + a(3)*b(0)

END SUBROUTINE QTNPRD

!*******************************************************************************!
!																				!
!	SUBROUTINE SYM_NEAR( isym, x, y, trmax)										!
!-------------------------------------------------------------------------------!
! 	Converts y(0:3) to the symmetry (set by isym) variant closest to x(0:3)		!
!...............................................................................!
!    symmetry codes:															!
!							0		triclinic									!
!							1		monoclinic									!
!							2		orthorhombic								!
!							3		trigonal									!
!							4		tetragonal									!
!							5		hexagonal									!
!							6		cubic										!
!																				!
!*******************************************************************************!
SUBROUTINE SYM_NEAR( isym, x, y, trmax)

	IMPLICIT NONE
	INTEGER(4)	isym, symon, i
	REAL(4)		x, y, a, b, asym, amin, trace, trmax
	DIMENSION	x(0:3), y(0:3), a(0:3), b(0:3), amin(0:3)

	REAL(4) QTNTCE

	COMMON/ SYMQCF/ asym(0:3, 1:30 ), symon(1:30,0:6)   ! symmetry operation data

	trmax= QTNTCE(x,y)						!	initialise cos(w/2)
	amin= y

	DO i= 2, 30
		IF ( symon(i, isym) .EQ. 1 ) THEN	!	check if symmetry in class

			a(0:3) = asym(0:3, i)
			CALL QTNPRD(a,y,b)

			IF ( b(0) .LT. 0.) b= -b		!	ensure in positive hemisphere
				
			trace= QTNTCE(x,b)				!	generate/test cos(w/2)

			IF ( trace .GT. trmax) THEN
				trmax= trace
				amin= b
			END IF
		END IF
	END DO

	y= amin

END SUBROUTINE SYM_NEAR

!*******************************************************************************!
!																				!
!   	BLOCK DATA SYMQTN														!
!-------------------------------------------------------------------------------!
!	Data for symmetry operations expressed as Euler symmetrics (0:3) and switch !
!	codes (1 for on, 0 for off) for their application in the 7 point groups.	!
!-------------------------------------------------------------------------------!
!    symmetry codes:															!
!							0		triclinic									!
!							1		monoclinic									!
!							2		orthorhombic								!
!							3		trigonal									!
!							4		tetragonal									!
!							5		hexagonal									!
!							6		cubic										!
!																				!
!*******************************************************************************!
BLOCK DATA SYMQTN

	IMPLICIT NONE
	INTEGER(4)	symon, i
	REAL(4)		asym

	COMMON/ SYMQCF/ asym(0:3, 1:30 ), symon(1:30,0:6)

	DATA (asym(i, 1),i= 0, 3)/ 1.0, 0.0, 0.0, 0.0 / 				!  identity

	DATA (asym(i, 2),i= 0, 3)/ 0.7071068, 0.7071068, 0.0, 0.0 / 	!  90 [ 1 0 0]
	DATA (asym(i, 3),i= 0, 3)/ 0.0,       1.0,   0.0, 0.0     / 	! 180 [ 1 0 0]
	DATA (asym(i, 4),i= 0, 3)/ 0.7071068,-0.7071068, 0.0, 0.0 / 	! 260 [ 1 0 0]

	DATA (asym(i, 5),i= 0, 3)/ 0.7071068, 0.0, 0.7071068, 0.0 / 	!  90 [ 0 1 0]
	DATA (asym(i, 6),i= 0, 3)/ 0.0,       0.0, 1.0, 0.0       / 	! 180 [ 0 1 0]
	DATA (asym(i, 7),i= 0, 3)/ 0.7071068, 0.0,-0.7071068, 0.0 / 	! 270 [ 0 1 0]

	DATA (asym(i, 8),i= 0, 3)/ 0.7071068, 0.0, 0.0, 0.7071068 / 	!  90 [ 0 0 1]
	DATA (asym(i, 9),i= 0, 3)/ 0.0,       0.0, 0.0, 1.0       / 	! 180 [ 0 0 1]
	DATA (asym(i,10),i= 0, 3)/ 0.7071068, 0.0, 0.0,-0.7071068 / 	! 270 [ 0 0 1]

	DATA (asym(i,11),i= 0, 3)/0.0, 0.7071068, 0.7071068, 0.0 / 		! 180 [ 1 1 0]
	DATA (asym(i,12),i= 0, 3)/0.0,-0.7071068, 0.7071068, 0.0 / 		! 180 [-1 1 0]

	DATA (asym(i,13),i= 0, 3)/0.0, 0.7071068, 0.0, 0.7071068 / 		! 180 [ 1 0 1]
	DATA (asym(i,14),i= 0, 3)/0.0,-0.7071068, 0.0, 0.7071068 / 		! 180 [-1 0 1]

	DATA (asym(i,15),i= 0, 3)/0.0, 0.0, 0.7071068, 0.7071068 / 		! 180 [ 0 1 1]
	DATA (asym(i,16),i= 0, 3)/0.0, 0.0,-0.7071068, 0.7071068 / 		! 180 [ 0-1 1]

!-----------------cubic {111} triads
	DATA (asym(i,17),i= 0, 3)/ 0.5, 0.5, 0.5, 0.5/					! 120 [ 1 1 1]
	DATA (asym(i,18),i= 0, 3)/ 0.5, -0.5,-0.5,-0.5/					! 240 [ 1 1 1]

	DATA (asym(i,19),i= 0, 3)/ 0.5,-0.5, 0.5, 0.5/					! 120 [-1 1 1]
	DATA (asym(i,20),i= 0, 3)/ 0.5, 0.5,-0.5,-0.5/					! 240 [-1 1 1]

	DATA (asym(i,21),i= 0, 3)/ 0.5, 0.5,-0.5, 0.5/					! 120 [ 1-1 1]
	DATA (asym(i,22),i= 0, 3)/ 0.5,-0.5, 0.5,-0.5/					! 240 [ 1-1 1]

	DATA (asym(i,23),i= 0, 3)/ 0.5, 0.5, 0.5,-0.5/					! 120 [ 1 1-1]
	DATA (asym(i,24),i= 0, 3)/ 0.5,-0.5,-0.5, 0.5/					! 240 [ 1 1-1]

!-------------hexgonal hexads, include 180 [ 0 0 1]: subset for trigonal
	DATA (asym(i,25),i= 0, 3)/ 0.866254, 0.0, 0.0, 0.5      / 		!  60 [ 0 0 1]
	DATA (asym(i,26),i= 0, 3)/ 0.5,      0.0, 0.0, 0.866254 / 		! 120 [ 0 0 1]
	DATA (asym(i,27),i= 0, 3)/ 0.5,      0.0, 0.0,-0.866254 / 		! 240 [ 0 0 1]
	DATA (asym(i,28),i= 0, 3)/ 0.866354, 0.0, 0.0,-0.5      / 		! 300 [ 0 0 1]

!--------------hexagonal diads, use with 180 [ 1 0 0]
	DATA (asym(i,29),i= 0, 3)/ 0.0,-0.5, 0.866254, 0.0      / 		! 180 [ 0 1-1 0]
	DATA (asym(i,30),i= 0, 3)/ 0.0,-0.5,-0.866254, 0.0      / 		! 180 [ 0-1 1 0]

!---------------------------------on/off codes for point groups
	DATA (symon(i,0),i=1,30)/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,	&
						&	 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0   /		!triclinic

	DATA (symon(i,1),i=1,30)/1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,	&
						&	 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0   /		!monoclinic

	DATA (symon(i,2),i=1,30)/1,0,1,0,0,1,0,0,1,0,0,0,0,0,0,	&
						&	 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0   /		!orthorhombic

	DATA (symon(i,3),i=1,30)/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,	&
						&	 0,0,0,0,0,0,0,0,0,0,1,1,0,0,0   /		!trigonal

	DATA (symon(i,4),i=1,30)/1,0,1,0,0,1,0,1,1,1,1,1,0,0,0,	&
						&	 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0   /		!tetragonal

	DATA (symon(i,5),i=1,30)/1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,	&
						&	 0,0,0,0,0,0,0,0,0,1,1,1,1,1,1   /		!hexagonal

	DATA (symon(i,6),i=1,30)/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,	&
						&	 1,1,1,1,1,1,1,1,1,0,0,0,0,0,0   /		!cubic

END BLOCK DATA SYMQTN
