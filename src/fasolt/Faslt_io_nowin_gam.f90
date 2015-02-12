!***********************************************************************!
!																		!
!		SUBROUTINE RDDATA												!
!-----------------------------------------------------------------------!
!		Reads data file for FASOLT3										!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE RDDATA(f_name)

	USE FAS_COM			! FASOLT main variables declaration
	USE FAS_CPL			! FASOLT material property block

	IMPLICIT NONE

	CHARACTER*(*) f_name
	LOGICAL(4)  l_file
	INTEGER(4)	i, j, k, i1, i2, i3, nde, npe

	INTEGER(4)	iomon, in, out, f_len
    COMMON/ iochan/ iomon, in, out
	in=52
	
!-------------------------common dialog to open file and get file name
	CALL INPFILE( in, l_file, f_name)		
	IF ( .NOT. l_file) CALL FAFERR(501)
!-----------------------------------------------------read FAFNER data

		READ(in) ninc	!	start increment
		ninc0=ninc+1
		READ(in) nlmnt, nnod, ngps, nsv	! size of mesh
		IF(		nlmnt.GT. maxels .OR. nnod.GT. maxnod .OR.	&
			&	ngps .GT. maxgps .OR. nsv .GT. maxsvs   ) CALL FAFERR(510)

!------------------------read element definitions
		DO i= 1, nlmnt
			READ(in) le(i), lt(i), nde, npe 
			READ(in)(ln(j,i), j= 1, nde)
			READ(in)(lg(j,i), j= 1, npe)
			nd(i)=nde
			np(i)=npe
		END DO

!----------------------read nodal variables 
		DO i= 1, nnod
			i3= 3*i
			i2= i3 - 1
			i1= i2 - 1
			READ(in) ib(i)
			READ(in) x(i1), x(i2), x(i3),  dx(i1), dx(i2), dx(i3)
			READ(in) fc(i1), fc(i2), fc(i3)
		END DO

!---------------------read integration point variables
		DO i= 1, ngps
			READ(in) eps(i), dep(i)
			READ(in) ( sv(j,i), j= 1, nsv)
			READ(in) (dsv(j,i), j= 1, nsv)
			READ(in) (ep(j,i), j= 1, 6)
			READ(in) (de(j,i), j= 1, 6)
			READ(in) (st(j,i), j= 1, 6)
			READ(in) (ds(j,i), j= 1, 6)
			READ(in) (rt(j,i), j= 1, 3)
			READ(in) (dr(j,i), j= 1, 3)
		END DO

!-----------------------------initialise material properies to zero
		e= 0.
		nss= 0
		ab= 0. ; an= 0.
		ampar= 0.
		smatrx= 0.

!------------------------------------------read material properties
		READ(in) nmatl, nmpar
		IF ( nmatl .GT. maxmds) CALL FAFERR(601)
		DO mat_code= 1, nmatl
			READ(in) nss(mat_code)
			IF ( nss(mat_code) .GT. maxss )		CALL FAFERR(602)
			IF ( nss(mat_code) .GT. maxsvs -3 )	CALL FAFERR(603)

!------------------------------------real elastic constants
			READ(in) (e(j, mat_code), j= 1, 9)

!-----------------------------------------read slip systems
			DO i= 1, nss(mat_code)
				READ(in) (ab(j,i,mat_code),j=1,3), &
					&  (an(j,i,mat_code),j=1,3)
			END DO
!---------------------------------------read ss properties
			DO i= 1, nss(mat_code)
				READ(in) (ampar(j,i,mat_code),j=1,nmpar)
			END DO

!--------------------------------read latent hardening matrix
			DO i=1,nss(mat_code)
				READ(in) (smatrx(j,i, mat_code),j=1,nss(mat_code))
			END DO

		END DO

!---------------------------read surface properties
		READ(in) amu, anu, aga

!---------------------------------read control data
		READ(in) ninc1, nxs, nfs
		DO i= 1, ninc1
			READ(in) tiem(i)
			IF ( nxs .GT. 0 ) THEN
				DO j= 1, nxs
					READ(in) ( x_spec(k,j,i), k= 1,15)
				END DO
			ENDIF
			IF ( nfs .GT. 0 ) THEN
				DO j= 1, nfs
					READ(in) ( f_spec(k,j,i), k= 1, 3)
				END DO
			ENDIF
		END DO

	CLOSE(in)

!----------------------------------------open report file
	OPEN( UNIT=iomon, FILE='fafrep.txt')


END SUBROUTINE RDDATA

!***********************************************************************!
!																		!
!		SUBROUTINE WTDATA												!
!-----------------------------------------------------------------------!
!		Writes results file for FASOLT3									!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE WTDATA(f_name)

	USE FAS_COM			! FASOLT main variables declaration
	USE FAS_CPL			! FASOLT material property block

	IMPLICIT NONE

	CHARACTER*(*) f_name
	LOGICAL(4)  l_file
	INTEGER(4)	i, j, k, i1, i2, i3, nde, npe, ierr


	INTEGER(4)	iomon, in, out
    COMMON/ iochan/ iomon, in, out

!----------------------------------------------------------open file
	OPEN( UNIT= out, FILE= f_name, FORM= 'UNFORMATTED',	&
			&	ACCESS= 'SEQUENTIAL', STATUS= 'UNKNOWN', IOSTAT= ierr)
			write(*,*) f_name
	IF ( ierr .NE. 0) CALL FAFERR(502)

!-----------------------------------------------write FAFNER results
		WRITE(out) inc	! increment

		WRITE(out) nlmnt, nnod, ngps, nsv	! size of mesh etc.

!------------------------write element definitions
		DO i= 1, nlmnt
			nde=nd(i)
			npe=np(i)
			WRITE(out) le(i), lt(i), nde, npe
			WRITE(out)(ln(j,i), j= 1, nde)
			WRITE(out)(lg(j,i), j= 1, npe)
		END DO

!----------------------write nodal variables 
		DO i= 1, nnod
			i3= 3*i
			i2= i3 - 1
			i1= i2 - 1
			WRITE(out) ib(i)
			WRITE(out) x(i1), x(i2), x(i3),  dx(i1), dx(i2), dx(i3)
			WRITE(out) fc(i1), fc(i2), fc(i3)
		END DO

!---------------------write integration point variables
		DO i= 1, ngps
			WRITE(out) eps(i), dep(i)
			WRITE(out) (sv(j,i), j= 1, nsv)
			WRITE(out) (dsv(j,i), j= 1, nsv)
			WRITE(out) (ep(j,i), j= 1, 6)
			WRITE(out) (de(j,i), j= 1, 6)
			WRITE(out) (st(j,i), j= 1, 6)
			WRITE(out) (ds(j,i), j= 1, 6)
			WRITE(out) (rt(j,i), j= 1, 3)
			WRITE(out) (dr(j,i), j= 1, 3)
		END DO

!------------------------------------------write material properties
		WRITE(out) nmatl, nmpar
		IF ( nmatl .GT. maxmds) CALL FAFERR(601)
		DO mat_code= 1, nmatl
			WRITE(out) nss(mat_code)

!------------------------------------write elastic constants
			WRITE(out) (e(j, mat_code), j= 1, 9)

!-----------------------------------------write slip systems
			DO i= 1, nss(mat_code)
				WRITE(out) (ab(j,i,mat_code),j=1,3), &
						&  (an(j,i,mat_code),j=1,3)
			END DO
!---------------------------------------write ss properties
			DO i= 1, nss(mat_code)
				WRITE(out) (ampar(j,i,mat_code),j=1,nmpar)
			END DO

!--------------------------------write latent hardening matrix
			DO i=1,nss(mat_code)
				WRITE(out) (smatrx(j,i, mat_code),j=1,nss(mat_code))
			END DO

		END DO

!---------------------------write surface properties
		WRITE(out) amu, anu, aga

!---------------------------------write control data
		WRITE(out) ninc1, nxs, nfs
		DO i= 1, ninc1
			WRITE(out) tiem(i)
			IF ( nxs .GT. 0 ) THEN
				DO j= 1, nxs
					WRITE(out) ( x_spec(k,j,i), k= 1,15)
				END DO
			ENDIF
			IF ( nfs .GT. 0 ) THEN
				DO j= 1, nfs
					WRITE(out) ( f_spec(k,j,i), k= 1, 3)
				END DO
			ENDIF
		END DO
!--------------------------write slip system activity
		WRITE(out) nss(mat_code)
		
		DO i=1,nss(mat_code)
			WRITE(out) gam(i)
		END DO



	CLOSE(out)

END SUBROUTINE WTDATA

!***********************************************************************!
!																		!
!		SUBROUTINE ITHEAD(ninc)											!
!-----------------------------------------------------------------------!
!		Writes iteration report header to monitor channel				!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE ITHEAD(inc)

	IMPLICIT NONE

	INTEGER(4)	iomon, in, out, inc

    COMMON/ iochan/ iomon, in, out

	WRITE(iomon,'(//,25X,''Increment number'',I4,/,10X,''Iteration'',2X,''Error(X)'',	&
					&	2X,''Error(F)'',2X,''Err(CPI)'',2X,''Decelrn.'',/)') inc
	WRITE(*    ,'(//,25X,''Increment number'',I4,/,10X,''Iteration'',2X,''Error(X)'',	&
					&	2X,''Error(F)'',2X,''Err(CPI)'',2X,''Decelrn.'',/)') inc

END SUBROUTINE ITHEAD

!***********************************************************************!
!																		!
!		SUBROUTINE ITEREP(nit,sdx,sfc,alp2)								!
!-----------------------------------------------------------------------!
!		Writes iteration status to monitor channel						!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE ITEREP(nit,sdx,sfc,cper,alp2)

	IMPLICIT NONE

	INTEGER(4)	iomon, in, out, nit
	REAL(4)		sdx, sfc, cper, alp2

    COMMON/ iochan/ iomon, in, out

	WRITE(iomon,'(12X,I3,4X,5(1X,E9.2))') nit, sdx, sfc, cper, alp2
	WRITE(*    ,'(12X,I3,4X,5(1X,E9.2))') nit, sdx, sfc, cper, alp2

END SUBROUTINE ITEREP

!***********************************************************************!
!																		!
!	SUBROUTINE POWER													!
!-----------------------------------------------------------------------!
!  Reports the power f.dx/dt for the the whole domain					!
!																		!
!***********************************************************************!
SUBROUTINE POWER

	USE FAS_COM

	IMPLICIT NONE
	INTEGER(4)	iomon, in, out

    COMMON/ iochan/ iomon, in, out

	WRITE(iomon,'(/,5X,''Power: '',G12.5)') DOT_PRODUCT( fc(1:3*nnod),dx(1:3*nnod))/tiem(ninc1)
	WRITE(*    ,'(/,5X,''Power: '',G12.5)') DOT_PRODUCT( fc(1:3*nnod),dx(1:3*nnod))/tiem(ninc1)


END SUBROUTINE

!***********************************************************************!
!																		!
!		SUBROUTINE INPFILE												!
!-----------------------------------------------------------------------!
!	Routine to OPEN FASOLT3 input file									!
!																		!
!-----------------------------------------------------------------------!
!Arguments:IN	INTEGER;	READ channel number 						!
!L_FILE		LOGICAL;    TRUE if file opened ok							!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!



SUBROUTINE 	INPFILE(in,l_file,f_name)
IMPLICIT 	none
LOGICAL(4) :: 		l_file
INTEGER ::			in,ierr
CHARACTER(100) ::	f_name,dir,dir1,dir2,dir3,f_name_old
CHARACTER(4) ::		ext

OPEN(in,FILE='dirname',STATUS='OLD', IOSTAT=ierr)
	IF ( ierr .EQ. 0 ) THEN
		READ(in,*) dir1,dir2,dir3
		dir=dir3
	ELSE
		WRITE(*,'("No dirname file found. Path is current directory.")')
		dir=''
	END IF
CLOSE(in)


WRITE(*,'(''Name of FASOLT input file (no extension): '',\)')
ext=".fsd"
CALL check_file(f_name,ext,dir)
OPEN(in ,FILE=TRIM(dir)//TRIM(f_name)//TRIM(ext) ,STATUS='OLD', FORM= 'UNFORMATTED',	&
		&	ACCESS= 'SEQUENTIAL', IOSTAT=ierr)
	IF ( ierr .EQ. 0 ) THEN
		l_file = .TRUE.
		ELSE
		WRITE(*,'("Could not open input file! ")')
		l_file = .FALSE.
	END IF
END SUBROUTINE INPFILE



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
CHARACTER(100) f_name_old
INTEGER :: i
LOGICAL :: here

DO i=1,3
READ(*,*) f_name
f_name_old=f_name
f_name=TRIM(dir)//TRIM(f_name)//TRIM(ext)
write(*,*) f_name
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
f_name=f_name_old
END SUBROUTINE check_file




