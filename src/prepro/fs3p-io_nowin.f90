
!***********************************************************************!
!																		!
!		SUBROUTINE WTDATA												!
!-----------------------------------------------------------------------!
!		Writes data file for FASOLT3									!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE WTDATA

	USE FAS_COM

	IMPLICIT NONE

	LOGICAL(4)  L_FILE
	INTEGER(4)	i, j, k, i1, i2, i3

	INTEGER(4), PARAMETER::	 IN= 51

!------------------------------------common dialog to open file
	CALL OUTFILE( IN, L_FILE)

!----------------------------------------write FAFNER data
		WRITE(in) ninc1	!	start increment

		WRITE(in) nlmnt, nnod, ngps, nsv	! size of mesh

!------------------------write element definitions
		DO i= 1, nlmnt
			WRITE(in) le(i), lt(i), nd(i), np(i)
			WRITE(in)(ln(j,i), j= 1, nd(i))
			WRITE(in)(lg(j,i), j= 1, np(i))
		END DO

!----------------------write nodal variables 
		DO i= 1, nnod
			i3= 3*i
			i2= i3 - 1
			i1= i2 - 1
			WRITE(in) ib(i)
			WRITE(in) x(i1), x(i2), x(i3),  dx(i1), dx(i2), dx(i3)
			WRITE(in) fc(i1), fc(i2), fc(i3)
		END DO

!---------------------write integration point variables
		DO i= 1, ngps
			WRITE(in) eps(i), dep(i)
			WRITE(in) ( sv(j,i), j= 1, nsv)
			WRITE(in) (dsv(j,i), j= 1, nsv)
			WRITE(in) (ep(j,i), j= 1, 6)
			WRITE(in) (de(j,i), j= 1, 6)
			WRITE(in) (st(j,i), j= 1, 6)
			WRITE(in) (ds(j,i), j= 1, 6)
			WRITE(in) (rt(j,i), j= 1, 3)
			WRITE(in) (dr(j,i), j= 1, 3)
		END DO

!------------------------------------------write material properties
		WRITE(in) nmatl, nmpar
		DO mat_code= 1, nmatl
			WRITE(in) nss(mat_code)

!------------------------------------write twinning limit
			WRITE(in) twinlim

!------------------------------------write twinning probability
			WRITE(in) twinprob

!------------------------------------write elastic constants
			WRITE(in) (e(j, mat_code), j= 1, 9)

!------------------------------------write slip system reorientation angles
			WRITE(in) (ssang(j, mat_code), j= 1, nss(mat_code))

!-----------------------------------------write slip systems
			DO i= 1, nss(mat_code)
				WRITE(in) (ab(j,i,mat_code),j=1,3), &
						&  (an(j,i,mat_code),j=1,3)
			END DO
!---------------------------------------write ss properties
			DO i= 1, nss(mat_code)
				WRITE(in) (ampar(j,i,mat_code),j=1,nmpar)
			END DO

!--------------------------------write latent hardening matrix
			DO i=1,nss(mat_code)
				WRITE(in) (smatrx(j,i, mat_code),j=1,nss(mat_code))
			END DO
		END DO
!---------------------------write surface properties
		amu=0.
		anu=0.
		aga=0.

		WRITE(in) amu, anu, aga


!---------------------------------write control data
		WRITE(in) ninc1, nxs, nfs
		DO i= 1, ninc1
			WRITE(in) tiem(i)
			IF ( nxs .GT. 0 ) THEN
				DO j= 1, nxs
					WRITE(in) ( x_spec(k,j,i), k= 1, 15)
				END DO
			ENDIF
			IF ( nfs .GT. 0 ) THEN
				DO j= 1, nfs
					WRITE(in) ( f_spec(k,j,i), k= 1, 3)
				END DO
			ENDIF
		END DO

	CLOSE(IN)

END SUBROUTINE WTDATA


!***********************************************************************!
!																		!
!		SUBROUTINE OUTFILE												!
!-----------------------------------------------------------------------!
!	Routine to OPEN file for writing        							!
!																		!
!-----------------------------------------------------------------------!
!	Arguments:	IN			INTEGER;	READ channel number 			!
!				L_FILE		LOGICAL;    TRUE if file opened ok			!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!

SUBROUTINE OUTFILE(IN,L_FILE)

IMPLICIT none
LOGICAL(4) ::  		l_file
INTEGER :: 			in, ierr
CHARACTER(100) ::	f_name,dir,dir1,dir2,dir3
CHARACTER(4) ::		ext

OPEN(in,FILE='dirname',STATUS='OLD', IOSTAT=ierr)
	IF ( ierr .EQ. 0 ) THEN
		READ(in,*) dir1,dir2,dir3
		dir=dir3
		WRITE(*,'(''FASOLT file directory: '',\)')
		WRITE(*,*) dir
	ELSE
		WRITE(*,'("No dirname file found. Path is current directory.")')
		dir='.'
	END IF
CLOSE(in)


WRITE(*,'(''Name of FASOLT input file (no extension): '',\)')
READ(*,*) f_name
ext=".fsd"
f_name=TRIM(dir)//TRIM(f_name)//TRIM(ext)

OPEN( UNIT= IN, FILE= F_NAME, FORM= 'UNFORMATTED',		&
	&	ACCESS= 'SEQUENTIAL', STATUS= 'UNKNOWN',IOSTAT= IERR)
		
IF ( ierr .EQ. 0 ) THEN
	l_file = .TRUE.
	ELSE
	WRITE(*,'("Could not open input file! ")')
	l_file = .FALSE.
END IF

END SUBROUTINE OUTFILE

!***********************************************************************!
!																		!
!		SUBROUTINE ORIFILE												!
!-----------------------------------------------------------------------!
!	Routine to OPEN file for reading ORI file         							!
!																		!
!-----------------------------------------------------------------------!
!	Arguments:	IN			INTEGER;	READ channel number 			!
!				L_FILE		LOGICAL;    TRUE if file opened ok			!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!

SUBROUTINE ORIFILE(in,l_file)

IMPLICIT none
LOGICAL(4) ::  		l_file
INTEGER :: 			in, ierr
CHARACTER(100) ::	f_name,dir,dir1,dir2,dir3
CHARACTER(4) ::		ext

OPEN(in,FILE='dirname',STATUS='OLD', IOSTAT=ierr)
	IF ( ierr .EQ. 0 ) THEN
		READ(in,*) dir1,dir2,dir3
		dir=dir2
		WRITE(*,'(''Orientation file directory: '',\)')
		WRITE(*,*) dir
	ELSE
		WRITE(*,'("No dirname file found. Path is current directory.")')
		dir='.'
	END IF
CLOSE(in)


WRITE (*,'(''Name of orientation data file (no extension): '',\)')
ext=".ori"
CALL check_file(f_name,ext,dir)

OPEN( UNIT= IN, FILE= F_NAME, FORM= 'FORMATTED',		&
			&	ACCESS= 'SEQUENTIAL', STATUS= 'OLD',IOSTAT= IERR)

IF ( ierr .EQ. 0 ) THEN
	l_file = .TRUE.
	ELSE
	WRITE(*,'("Could not open orientations file! ")')
	l_file = .FALSE.
END IF

END SUBROUTINE ORIFILE

!***********************************************************************!
!																		!
!		SUBROUTINE MATFILE												!
!-----------------------------------------------------------------------!
!	Routine to OPEN file for reading material properties file         							!
!																		!
!-----------------------------------------------------------------------!
!	Arguments:	IN			INTEGER;	READ channel number 			!
!				L_FILE		LOGICAL;    TRUE if file opened ok			!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!

SUBROUTINE MATFILE(in,l_file)
IMPLICIT none
LOGICAL(4) ::  		l_file
INTEGER :: 			in, ierr
CHARACTER(100) ::	f_name,dir,dir1,dir2,dir3
CHARACTER(4) ::		ext


OPEN(in,FILE='dirname',STATUS='OLD', IOSTAT=ierr)
	IF ( ierr .EQ. 0 ) THEN
		READ(in,*) dir1,dir2,dir3
		dir=dir1
		WRITE(*,'(''Material data file directory: '',\)')
		WRITE(*,*) dir
	ELSE
		WRITE(*,'("No dirname file found. Path is current directory.")')
		dir='.'
	END IF
CLOSE(in)


WRITE (*,'(''Name of material data file (no extension): '',\)')
ext=".fmt"

CALL check_file(f_name,ext,dir)

OPEN( UNIT= IN, FILE= F_NAME, FORM= 'FORMATTED',		&
			&	ACCESS= 'SEQUENTIAL', STATUS= 'OLD',IOSTAT= IERR)

IF ( ierr .EQ. 0 ) THEN
	l_file = .TRUE.
	ELSE
	WRITE(*,'("Could not open materials file! ")')
	l_file = .FALSE.
END IF
		
END SUBROUTINE MATFILE

!***********************************************************************!
!																		!
!		SUBROUTINE check_file											!
!-----------------------------------------------------------------------!
!	Routine to help with inputing file to open							!
!																		!
!-----------------------------------------------------------------------!
!Arguments:f_name ext dir; to define file and path						!
!here		LOGICAL;    TRUE if file exhists							!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!

SUBROUTINE check_file(f_name,ext,dir)
IMPLICIT none
CHARACTER(100), INTENT(OUT) :: f_name
CHARACTER(100), INTENT(IN) :: dir
CHARACTER(4), INTENT(IN) :: ext
INTEGER :: i
LOGICAL :: here

DO i=1,3
READ(*,*) f_name
f_name=TRIM(dir)//TRIM(f_name)//TRIM(ext)

	INQUIRE(FILE=f_name, EXIST=here)

	IF (here) THEN
		EXIT
	ELSE IF (i .LT. 3) THEN
		WRITE (*,'("File not found. Try again.")')
		CYCLE
	ELSE 
		WRITE (*,'("File not found. Exiting.")')
	STOP
	END IF
END DO

END SUBROUTINE check_file
