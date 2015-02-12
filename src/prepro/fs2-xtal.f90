!*******************************************************************************!
!																				!
!		SUBROUTINE PERTID														!
!        Perturbation of ideal orientation										!
!																				!
!*******************************************************************************!
SUBROUTINE PERTID( ph1, phi, ph2, chi, r1, r2, r3 )

!	USE DFLIB

	IMPLICIT NONE
	REAL(4)		ph1, phi, ph2, chi, r1, r2, r3, a, b
	DIMENSION	a(1:3, 1:3), b( 1:3, 1:3)

!----------------------------------------ideal and perturbation matrices
	CALL EUATOA( ph1, phi, ph2, a) 

	CALL RANDOM_NUMBER(r1)
	CALL RANDOM_NUMBER(r2)
	CALL RANDOM_NUMBER(r3)
	r1= ACOS( 2.*r1 -1.)
	r2= 6.28318530718 *r2
	r3= chi*r3
	CALL AXATOA( r3, r1, r2, b)

!--------------------------------------------------multiply, get angles
	a= MATMUL(a,b)


	IF ( ABS (a(3,3)) .GE. 1.) then
		r1= ATAN2( a(1,2), a(1,1) )/2.
		r2= 0.
		r3= r1
	ELSE
		r1= ATAN2( a(3,1), -a(3,2) )
		r2= ACOS(  a(3,3) )
		r3= ATAN2( a(1,3),  a(2,3) )
	END IF

END SUBROUTINE

!***************************************************************************!
!																			!
!         SUBROUTINE EUATOA(ph1,phi,ph2,a) 									!
!---------------------------------------------------------------------------!
!	  Orientation matrix, a, from Euler angles ph1, phi, ph2(in radians)	!
!     with Bunge's definition ( all anticlockwise, about current z, x, z).  !
!	     Matrix has crystal definition of specimen axes as columns.			!
!																			!
!***************************************************************************!
SUBROUTINE EUATOA( ph1, phi, ph2, a)

	IMPLICIT NONE

	REAL(4)		ph1, phi, ph2, a
	REAL(4)		c, s, c1, c2,s1, s2, c1c2, c1s2, s1c2, s1s2

	DIMENSION	a(1:3,1:3)

!-----------------------------Form trig function terms
	c1= COS(ph1)
	s1= SIN(ph1)
	c2= COS(ph2)
	s2= SIN(ph2)
	c=  COS(phi)
	s=  SIN(phi)
	c1c2=c1*c2
	s1s2=s1*s2
	c1s2=c1*s2
	s1c2=s1*c2

!---------------------------Form the transform matrix
	a(1,1)= c1c2 - s1s2*c
	a(2,1)=-c1s2 - s1c2*c
	a(3,1)= s1*s
	a(1,2)= s1c2 + c1s2*c
	a(2,2)=-s1s2 + c1c2*c
	a(3,2)=-c1*s
	a(1,3)= s2*s
	a(2,3)= c2*s
	a(3,3)= c


END SUBROUTINE EUATOA

!***************************************************************************!
!																			!
!      SUBROUTINE AXATOA( omega, theta, psi, a)								!
!---------------------------------------------------------------------------!
!	Routine to calculate the orientation matrix, a, from axis/angle	measure !
!	with angle omega, about the axis with polar angle theta and azimuthal	!
!	angle psi. All angles in radians.										!
!	     Matrix has crystal definition of specimen axes as columns.			!
!																			!
!***************************************************************************!
SUBROUTINE AXATOA(omega,theta,psi,a)

	IMPLICIT NONE

	INTEGER(4)	i, j, k 
	REAL(4)		omega, theta, psi, a, d, st, cw, sw	

	DIMENSION	a(1:3,1:3), d(3)

!-------------------------------------------Rotation axis
    st= SIN(theta)
    d(1)= COS(psi)*st
    d(2)= SIN(psi)*st
    d(3)= COS(theta)

!-------------------------------Form matrix in three stages
    cw=COS(omega)
    sw=SIN(omega)

	a= 0.

    DO  i= 1, 3
		a(i,i)= cw
		DO  J= 1, 3
			a(i,j)= a(i,j) + (1.-cw)*d(i)*d(j)
		END DO

		j= 1 + MOD(i,3)
		k= 1 + MOD(j,3)
		a(j,k)=a(j,k)-sw*d(i)
		a(k,j)=a(k,j)+sw*d(i)

	END DO

END SUBROUTINE AXATOA
