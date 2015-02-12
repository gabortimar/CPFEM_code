!*******************************************************************************!
!																				!
!	SUBROUTINE SETSTATE															!
!		Sets the internal element state											!
!																				!
!*******************************************************************************!
SUBROUTINE SETSTATE( nein, location)

	USE FAS_COM
!	USE DFLIB

	IMPLICIT NONE
	LOGICAL(4)			l_file
	INTEGER(4)	nein, iselect, i, j, k, location, ilocal, set_type, nclusters, iclus
	REAL(4)		ph1, phi, ph2, ch1, chi, ch2,  r1, r2, r3, dummy
	DIMENSION	nein(1:3), location(1:3, 1:*), ilocal(1:3)

	REAL(4)		sumf, r
	REAL(4), DIMENSION(1:12)::		frac


	INTEGER(4), ALLOCATABLE:: cluster(:)
	REAL(4),	ALLOCATABLE:: ass_ang(:,:), ass_gps(:,:)

	ALLOCATE ( cluster(1:nlmnt), ass_ang(1:3, 1:nlmnt), ass_gps(1:3, 1:ngps) )

!-------------------------------------get number of slip systems
	WRITE(*,'(//,2X,''Enter max. no. of slip systems involved: '',\)') 
	READ(*,*) nsv
	nsv= nsv + 3

!------------------------------------------get type of operation
	WRITE(*,'(//,   ''  Select bicrystal         (1)'',/,		&
				&	''         perturbed-ideal   (2)'',/,	&
				&	''         included grain    (3)'',/,	&
				&	''         file orientations (4)'',/,	&
				&	''               -or random  (5): '',\)') 
	READ(*,*) iselect

	SELECT CASE( iselect)

		CASE(1)
			WRITE(*,'(10X,''Enter Euler angles for bottom half: '',\)')
			READ(*,*) ph1,phi,ph2
			ph1= ph1/ 57.295779513
			phi= phi/ 57.295779513
			ph2= ph2/ 57.295779513

			WRITE(*,'(10X,''Enter Euler angles for top half: '',\)')
			READ(*,*) ch1,chi,ch2
			ch1= ch1/ 57.295779513
			chi= chi/ 57.295779513
			ch2= ch2/ 57.295779513

		CASE(2)
			WRITE(*,'(10X,''Enter Euler angles and chi spread: '',\)')
			READ(*,*) ph1,phi,ph2,chi
			ph1= ph1/ 57.295779513
			phi= phi/ 57.295779513
			ph2= ph2/ 57.295779513
			chi= chi/ 57.295779513

		CASE(3)
			WRITE(*,'(10X,''Enter Euler angles for central region: '',\)')
			READ(*,*) ph1,phi,ph2
			ph1= ph1/ 57.295779513
			phi= phi/ 57.295779513
			ph2= ph2/ 57.295779513

		CASE(4)
			WRITE(*,'(10X,''Enter number of orientation sets: '',\)')
			READ(*,*) nmatl

		CASE DEFAULT
				! no i/o required

	END SELECT


!	CALL SEED( RND$TIMESEED )		!	initialise random sequence !DIGITAL Fortran only?

!--------get list of different location orientations ( for options 2, 4 & default)


!	IF ( iselect.NE. 1 .AND. iselect.NE. 3 .AND. iselect.NE. 4 ) THEN
	IF ( iselect.NE. 1 .AND. iselect.NE. 3 ) THEN

		ass_ang(1,1:nlmnt)= -999.	!	flag for not set
		set_type= -1				!	modular number for element code

		DO i= 1, nlmnt

			IF ( ass_ang(1,i) .EQ. -999.) THEN		!	new 'grain'

				set_type=  set_type + 1		! increment type number

!-------------------------------get next orientation depending on option
				SELECT CASE( iselect)

					CASE(2)				!	perturbed ideal  orientations
						CALL PERTID( ph1, phi, ph2, chi, r1, r2, r3 )


					CASE DEFAULT		!	random orientations
						CALL RANDOM_NUMBER(r1)
						CALL RANDOM_NUMBER(r2)
						CALL RANDOM_NUMBER(r3)
						r1= 6.2832*r1
						r2= ACOS(2.*r2-1.)
						r3= 6.2832*r3

				END SELECT

!-----------check here and upstream for same location, set angles and code
				ilocal(1:3)=location(1:3,i)	! location indices for i-th element

				DO j= i, nlmnt
					IF (  ( ilocal(1) .EQ. location(1,j)) .AND.	&
						& ( ilocal(2) .EQ. location(2,j)) .AND.	&
						& ( ilocal(3) .EQ. location(3,j)) ) THEN
							ass_ang(1,j)= r1
							ass_ang(2,j)= r2
							ass_ang(3,j)= r3

							le(j)= 1 + set_type
					END IF
				END DO


			END IF
		END DO
				
	END IF


!------------------------------------------------------Loop on elements setting orientation
	IF ( iselect .EQ. 4 ) THEN		! file orientations

		WRITE(*,'(/,15X,''Assigning orientation sets: enter fraction for: '')')
		DO i= 1, nmatl
			WRITE(*,'(18X,''set '',I3,'': '',\)') i
			READ(*,*) frac(i)
		END DO

!------------------------------normalise and make cumulative
		sumf= SUM(frac(1:nmatl)) ; frac(1:nmatl)= frac(1:nmatl)/sumf
		IF ( nmatl.GT.1) THEN
			DO i=2, nmatl
				frac(i)= frac(i)+frac(i-1)
			END DO
		END IF

!-----------------------------assign set to clusters (random)
		cluster= le		! assign element cluster number

		nclusters=0
		DO i= 1, nlmnt
			IF ( cluster(i) .GT. nclusters) nclusters= cluster(i)
		END DO

		DO i= 1, nclusters
			CALL RANDOM_NUMBER(r)
			k=0
			DO j= 1, nmatl
				IF ( r .GT. frac(j) ) k=j
			END DO

			DO j= 1, nlmnt
				IF ( cluster(j) .EQ. i ) le(j)= k+1
			END DO
		END DO


	!	DO i= 1, nlmnt
	!		CALL RANDOM_NUMBER(r)
	!		k=0
	!		DO j= 1, nmatl
	!			IF ( r .GT. frac(j) ) k=j
	!		END DO
	!		le(i)= k+1	
	!	END DO

		DO j= 1, nmatl

			WRITE(*,'(25X,''<enter> to input orientation for set'',I3,\)') j
	!		READ(*,*) 
	!		CALL ORIFILE(51,l_file)

	!			DO i=1, nlmnt
	!				IF ( le(i) .EQ. j ) THEN
	!					IF ( cluster(i) .GT. 0 ) THEN

	!						READ(51,*) dummy, r1, r2, r3

	!						iclus= cluster(i)

	!						DO k= i, nlmnt
	!							IF ( cluster(k) .EQ. iclus ) THEN
	!								ass_ang(1,k)= r1	! assign Euler angles to k-th element
	!								ass_ang(2,k)= r2
	!								ass_ang(3,k)= r3
	!								cluster(k)= -cluster(k)
	!							END IF
	!						END DO

	!					END IF
	!				END IF
	!			END DO
		
			!	DO i=1, nlmnt
			!		IF ( le(i) .EQ. j ) THEN
			!			READ(51,*) dummy, r1, r2, r3
			!			ass_ang(1,i)= r1
			!			ass_ang(2,i)= r2
			!			ass_ang(3,i)= r3
			!		END IF
			!	END DO

			CLOSE(51)

			! Read orientations again, this time for every GP:---------------------------
			READ(*,*) 	! just a quick fix, to help read ori file twice.
			CALL ORIFILE(51,l_file)
			do i = 1, ngps
				READ(51,*) dummy, r1, r2, r3
				ass_gps(1,i)= r1
				ass_gps(2,i)= r2
				ass_gps(3,i)= r3
			end do
			CLOSE(51)
			! GP-wise orientations read.-------------------------------------------------


		END DO

	ELSE

		DO i=1,nlmnt

			SELECT CASE( iselect)

				CASE(1)		!	allocate bottom and top halves, set types
					IF ( location(3,i) .LE. nein(3)/2 ) THEN
						r1=ph1
						r2=phi
						r3=ph2
						le(i)= 1

					ELSE
						r1=ch1
						r2=chi
						r3=ch2
						le(i)= 2
					END IF

					ass_ang(1,i)= r1
					ass_ang(2,i)= r2
					ass_ang(3,i)= r3

				CASE(3)		!  'included' grains
					IF (	( ABS( 2*location(3,i) -nein(3) -1 ).LT. nein(3)/2 ) .AND.	&
						&	( ABS( 2*location(2,i) -nein(2) -1 ).LT. nein(2)/2 ) .AND.	&
						&	( ABS( 2*location(1,i) -nein(1) -1 ).LT. nein(1)/2 ) )  THEN
						r1=ph1
						r2=phi
						r3=ph2
						le(i)= 1

					ELSE
						CALL RANDOM_NUMBER(r1)
						CALL RANDOM_NUMBER(r2)
						CALL RANDOM_NUMBER(r3)
						r1=6.2832*r1
						r2=acos(2.*r2-1.)
						r3=6.2832*r3
						le(i)= 2

					END IF
				
					ass_ang(1,i)= r1
					ass_ang(2,i)= r2
					ass_ang(3,i)= r3

			END SELECT
		END DO
	END IF

!------------------------------------------assign orientation state variable values
	DO i= 1, nlmnt
		DO j=1,np(i)
			k= np(i)*(i-1) + j
			!sv(1,k)=ass_ang(1,i)	! assign same Euler angles to every GP inside element
			!sv(2,k)=ass_ang(2,i)
			!sv(3,k)=ass_ang(3,i)
			sv(1,k)=ass_gps(1,k)	! assign Euler angles from GP-wise list
			sv(2,k)=ass_gps(2,k)
			sv(3,k)=ass_gps(3,k)
		END DO

	END DO

!---------------------------------------------------------null for all other states
	DO i=1,ngps
		eps(i)=0.
		dep(i)=0.
		st(1:6,i)=0.
		ep(1:6,i)=0.
		ds(1:6,i)=0.
		de(1:6,i)=0.
		rt(1:3,i)=0.
		dr(1:3,i)=0.
	END DO


	DEALLOCATE ( ass_ang, cluster)


END SUBROUTINE
