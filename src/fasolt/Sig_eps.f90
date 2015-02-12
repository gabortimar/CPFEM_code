!***********************************************************************!
!																		!
!   sig_eps: Routines for strain and rotations, stress and ds/de 		!
!																		!
!***********************************************************************!

!***********************************************************************!
!																		!
!        SUBROUTINE DSTRNS(np,sfd,dx,de,r,deff)							!
!-----------------------------------------------------------------------!
!		- Linear approximations to the strains and plane rigid 			!
!		  body rotation increments at an integration point.				!
!-----------------------------------------------------------------------!
!  Arguements   npl    INTEGER;			No. of s.f.d's (nodes) involved	!
!               sfd   REAL ARRAY(1:3,1:NP);		S.F.D's					!
!               dxl    REAL ARRAY(1:NP,1:3);	Displacements at nodes  !
!               del    REAL ARRAY(1:6);			Strain increments       !
!               rl     REAL ARRAY(1:3);			Rotation increments		!
!               deff  REAL;		        von Mises effective value of de	!
!-----------------------------------------------------------------------!
!  Note conventions:  r; anticlockwise about X, Y, Z					!
!					  de;  11, 22, 33, 23, 31, 12	Shears doubled		!
!																		!
!***********************************************************************!
SUBROUTINE DSTRNS(npl, sfd, dxl, del, rl, deff)

	USE FAS_COM			! FASOLT main variables declaration

	IMPLICIT NONE

	INTEGER(4)	npl
	REAL(4)		sfd, dxl, del, rl, deff, ddxl
	INTEGER(4)	i, j, k
	REAL(4)		sum, FJ2VOM

	DIMENSION	sfd( 1:3, 1:npl), dxl( 1:3, 1:npl), del( 1:6), rl( 1:3)
	DIMENSION   ddxl( 1:3, 1:3)


    ddxl= MATMUL( sfd, TRANSPOSE(dxl))	!	Displacement gradients

	sum= 0.
	DO i= 1, 3
		j= 1 + mod(i,3)
		k= 1 + mod(j,3)
!------------------------------------------strains
		del(i)= ddxl(i,i)
		del(i+3)= ( ddxl(k,j) + ddxl(j,k) )
!------------------------------------------rotations
		rl(i)=    ( ddxl(k,j) - ddxl(j,k) )/2.
		sum= sum +rl(i)*rl(i)
	END DO

!----------von Mises effective value of increment and test
    deff= FJ2VOM(.FALSE.,del)

!----------------------------------------test magnitudes
	IF(  sum .GT. smalr) CALL FAFERR(905)
	IF( deff .GT. smale) CALL FAFERR(906)

END SUBROUTINE DSTRNS

!***********************************************************************!
!																		!
!      SUBROUTINE DSTRSS(iel,jgp,deff,es,de,st,ds,t,r,d)				!
!-----------------------------------------------------------------------!
!	-Determination of stresses and optional modulii along with			!
!	rotational modifications to stress and strain. Stresses are			!
!	calculated for crystal plasticity, the modulii are derivatives		!
!	for use in Newton-Raphson scheme.									!
!-----------------------------------------------------------------------!
!	Arguments:	iel		INTEGER,		Element property code			!
!               jgp		INTEGER,		Integration point number		!
!               deff    REAL,		Increment effective plastic strain	!
!               es,de   REAL ARRAYS(1:6),Total strains, strain incs.,	!
!               st,ds					 total stresses and stress		!
!										 incs. Order: 11,22,33,23,31,12	!
!               r       REAL ARRAY(1:3), R.B. Rotation (anticl. 1,2,3)  !
!               chi     REAL ARRAY(1:6, 1:6), Modulii					! 
!               t       REAL ARRAY(1:6), At return temp. updated stress !
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE DSTRSS( iel, jgp, deff, es, de, st, ds, t, r, chi)

	USE FAS_CPL		!material property module

	IMPLICIT NONE

	INTEGER(4)	iel, jgp 
	REAL(4)		deff, d

	INTEGER(4)	i, j, ifail
	REAL(4)		sinit, einit
	REAL(4)		FJ2VOM

    	REAL(4), DIMENSION( 1:3)::		r
	REAL(4), DIMENSION( 1:6)::		es, de, st, ds, t, ep, s
	REAL(4), DIMENSION( 1:6, 1:6):: 	chi, el_comp
	
	mat_code= ABS(iel)		!  material property set pointer

	CALL ELASTC( jgp, el_comp)	!	Get elastic compliances

	CALL TENS_ROT(r,st,t)		!	Total start stresses in updated frame

    s= t + ds	!	Initial guess plastic state using existing stress
    sinit= FJ2VOM(.TRUE., s)
    einit= FJ2VOM(.FALSE.,de)


!-----------------------------------------Test for elastic only solution
	IF ( iel.GT.0 .AND. sinit.GT.0. .AND. einit.GT.0. ) THEN	!Plastic

		CALL CPSTRS(jgp,el_comp,de,r,s,t,ep,chi)

!--------Transfer s to t stress,get eff. strain increment
		t= s
		deff=FJ2VOM(.FALSE.,ep)

!--------------Symmetrised elasto-plastic matrix
		DO i= 1, 6
			DO j= i, 6
				d=( chi(i,j) + chi(j,i) )/2.
				chi(i,j)=d
				chi(j,i)=d
			END DO
		END DO

	ELSE											!	Elastic only solution

		deff=0.

		CALL FACTOR(6,el_comp,ifail)
		DO i=1,6
			ep(1:6)=0.
			ep(i) = 1.
			CALL SOLVER(6,el_comp,ep)
			DO j=1,6
				chi(j,i)= ep(j)
			END DO
		END DO

		t= MATMUL( chi, de)

	END IF
	
	ds= t - st		!	Stress increments


END SUBROUTINE DSTRSS

!***********************************************************************!
!																		!
!		SUBROUTINE TENS_ROT(r,a,b)										!
!-----------------------------------------------------------------------!
!	Routine to rotate symmetrical 3x3 tensor, in compact notation.		!
!-----------------------------------------------------------------------!
!	Arguments:		r  REAL ARRAY(1:3); small rotations about x,y,z		!
!					a  REAL ARRAY(1:6); initial basis tensor			!
!					b  REAL ARRAY(1:6); rotated basis tensor			!
!-----------------------------------------------------------------------!
!    NB  The condensed tensor is tensor form i.e. not doubled shears	!
!																		!
!***********************************************************************!
SUBROUTINE TENS_ROT( r, a, b)

	IMPLICIT NONE

	INTEGER(4)  i, j, k
	REAL(4)		r, a, b,  rot, old, new

	DIMENSION	r( 1:3), a( 1:6), b( 1:6)
	DIMENSION   	rot( 1:3, 1:3), old( 1:3, 1:3), new( 1:3, 1:3) 

!------------------------------------allocate to full arrays
	DO i= 1, 3
		j= 1 + MOD(i,3)
		k= 1 + MOD(j,3)

		old(i,i)= a(i)
		old(j,k)= a(i+3)
		old(k,j)= old(j,k)

		rot(j,k)=  r(i)
		rot(k,j)= -r(i)
		rot(i,i)= SQRT( 1. - r(j)*r(j) -r(k)*r(k) )
	END DO

!---------------------------------------rotate basis
	new= MATMUL(TRANSPOSE(rot), MATMUL(old, rot) )

!------------------------------------condense result
	DO i= 1, 3
		j= 1 + MOD(i,3)
		k= 1 + MOD(j,3)

		b(i)  = new(i,i)
		b(i+3)= new(j,k)
	END DO

END SUBROUTINE TENS_ROT

!***********************************************************************!
!																		!
!        REAL FUNCTION FJ2VOM( lstress, x)								!
!																		!
!-----------------------------------------------------------------------!
!       Von Mises effective value, for variables which can include		!
!       hydrostatic, of either stress or strain (eng. shear strain)		!
!-----------------------------------------------------------------------!
!    Arguments	lstress	LOGICAL;		TRUE if stresses				!
!               x		REAL ARRAY(1:6); Stress or strain values		!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
REAL FUNCTION FJ2VOM( lstress, x)

	IMPLICIT NONE

	LOGICAL(4)	lstress
	INTEGER(4)	i, j
	REAL(4)		x, sum, amult

    DIMENSION	x( 1:6)

    IF ( lstress ) THEN
         amult= 6.0
    ELSE
         amult= 1.5
    ENDIF

    sum=0.
    DO i= 1, 3
		j= 1 + MOD(i,3)
        sum= sum +(x(i)-x(j))*(x(i)-x(j)) + amult*x(i+3)*x(i+3)
	END DO

    FJ2VOM= SQRT( sum /4.5)

END FUNCTION FJ2VOM
