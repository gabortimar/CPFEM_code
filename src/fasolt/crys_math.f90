!***********************************************************************!
!																		!
!                       Fasolt3   -crys_math							!
!                        Pete Bate 2000									!
!-----------------------------------------------------------------------!
!                  Crystal routines for 3-d FEM				 			!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
!
!!	SLPTAU(SIG,B,TAU) REPLACABLE BY  TAU= MATMUL(TRANSPOSE(B),SIG)

!!	SLPEPS(EPS,B,GAM) REPLACABLE BY  EPS= MATMUL(B,GAM)
 
!
!***********************************************************************!
!																		!
!         SUBROUTINE CTMATRX(ph1,phi,ph2,a) 							!
!-----------------------------------------------------------------------!
!		Orientation matrix, a, from Euler angles (in radians) and		!
!		augment terms used for angle change calculations.				!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE CTMATRX( ph1, phi, ph2, a, aug)

	IMPLICIT NONE

	REAL(4)		ph1, phi, ph2, a, aug
	REAL(4)		c, s, c1, c2,s1, s2, c1c2, c1s2, s1c2, s1s2

	DIMENSION	a(1:3,1:3), aug( 1:3)

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

!------Extra terms for orientation change calculations
	aug(1)=s2
	aug(2)=c2
	aug(3)=s

END SUBROUTINE CTMATRX

!***********************************************************************!
!																		!
!      SUBROUTINE BTMATRX( a)											!
!-----------------------------------------------------------------------!
!  Routine for the slip rate to strain rate and spin matrices			!
!   in the current orientation, given by three Euler angles.			!
!-----------------------------------------------------------------------!
!  Note conventions:  r; anticlockwise about X, Y, Z					!
!					  de;  11, 22, 33, 23, 31, 12						!
!																		!
!***********************************************************************!
SUBROUTINE BTMATRX( a )

	USE FAS_CPL

	IMPLICIT NONE

	REAL(4)		a
	INTEGER(4)	i, j, l
	REAL(4)		sdr, spn, sum1, sum2

    	DIMENSION	a( 1:3, 1:3), sdr( 1:3), spn( 1:3)

!------------------------------Form the slip rate to strain rate matrix
	DO i=1,nss(mat_code)
		DO  j=1,3
			sum1= 0.
			sum2= 0.
			DO  l=1,3
				sum1= sum1 + a(l,j)*ab(l,i,mat_code)
				sum2= sum2 + a(l,j)*an(l,i,mat_code)
			END DO
			sdr(j)=sum1
			spn(j)=sum2
		END DO

		b(1,i)= sdr(1)*spn(1) 
		b(2,i)= sdr(2)*spn(2)
		b(3,i)= sdr(3)*spn(3)
		b(4,i)= sdr(2)*spn(3) + sdr(3)*spn(2)
		b(5,i)= sdr(3)*spn(1) + sdr(1)*spn(3)
		b(6,i)= sdr(1)*spn(2) + sdr(2)*spn(1)

!------------------------------and the slip rate to global spin matrix

		rot(1,i)=( sdr(2)*spn(3)-sdr(3)*spn(2) )/2.
		rot(2,i)=( sdr(3)*spn(1)-sdr(1)*spn(3) )/2.
		rot(3,i)=( sdr(1)*spn(2)-sdr(2)*spn(1) )/2.

	END DO

END SUBROUTINE BTMATRX

!***********************************************************************!
!																		!
!       SUBROUTINE ECHROT( a, aug, r, dp1, dph, dp2)					!
!       P. Bate  1987													!
!-----------------------------------------------------------------------!
!     Routine to give the rate of change of Euler angles of an			! 
!    orientation as the result of a spin referred to global frame.		!
!-----------------------------------------------------------------------!
!   Arguments:															!
!			a			REAL ARRAY(3,3); Orientation matrix				!
!			aug			REAL ARRAY(1:3); Terms for Euler angle change	!
!           r			REAL ARRAY(1:3); Rigid body spins about x,y,z	!
!           dp1,dph,dp2 REALs; changes of Euler angles					!
!-----------------------------------------------------------------------!
!   Convention:  eps(1:6)  e(11),e(22),e(33),e(12),e(23),e(31)			!
!                rtn(1:3)    23,31,12 (anticlock. about x,y,z)			!
!																		!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE ECHROT( a, aug, r, dp1, dph, dp2)

	USE FAS_CPL

	IMPLICIT NONE

	REAL(4)		a, aug, r, dp1, dph, dp2
	REAL(4)		rtn,rx,small

	DIMENSION	a( 1:3, 1:3), aug( 1:3), r( 1:3)
    	DIMENSION	rtn( 1:3), rx( 1:3)


	DATA small/ 1.e-5/
  
	rtn= MATMUL( rot(1:3,1:nss(mat_code)), gam(1:nss(mat_code))) ! spin due to slip
    	rtn= rtn - r				! subtract deformation r-b spin
	rx= MATMUL( a, rtn)			! transform to crystal values

!---------------------------Euler angle rates with 'gimbal lock' to Ph1
	IF ( ABS( aug(3) ) .LT. small) THEN
		dp1=rtn(3)
	ELSE
		dp1=( aug(1)*rx(1) + aug(2)*rx(2) )/aug(3)
	ENDIF

	dph= aug(2)*rx(1) - aug(1)*rx(2)
	dp2= rx(3) - dp1*a(3,3)


END SUBROUTINE ECHROT

