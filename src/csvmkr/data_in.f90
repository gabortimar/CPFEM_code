!***********************************************************************!
!																		!
!		SUBROUTINE FPREAD												!
!-----------------------------------------------------------------------!
!		Reads data file for FASPLOT3									!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE FPREAD(f_name)

	USE FAS_COM

	IMPLICIT NONE

	LOGICAL(4) l_file
	INTEGER(4)  i, j, k, i1, i2, i3, nde, npe, res, twinlim
	REAL(4)		twinprob
	
	INTEGER(4), DIMENSION(1:20):: nss
	
    	CHARACTER(100) ::	f_name
	INTEGER(4), PARAMETER::  in= 51

!-----open file to read here-------!!!!!!!!


		CALL INPFILE( in, l_file)		
		!IF ( .NOT. l_file) CALL FAFERR(501) 

		READ(in, ERR=991) ninc	!	 increment number
		ninc0=ninc+1

		READ(in,ERR=991) nlmnt, nnod, ngps, nsv	! number of elements, nodes, gauss points; nsv: number of slip systems + 3
		WRITE (*,*) nlmnt,maxels,nnod,maxnod,ngps,maxgps
		IF( nlmnt.GT.maxels .OR. nnod.GT.maxnod .OR. ngps.GT.maxgps) STOP

!------------------------read element definitions
		DO i= 1, nlmnt
			READ(in,ERR=991) le(i), lt(i), nde, npe
			READ(in,ERR=991)(ln(j,i), j= 1, nde)
			READ(in,ERR=991)(lg(j,i), j= 1, npe)
			nd(i)=nde
			np(i)=npe
		END DO

!----------------------exclude zero integration point ( tie) elements
		i=1
		DO WHILE ( i .LE. nlmnt)
			IF (np(i).EQ.0) THEN
				nlmnt= nlmnt-1
				IF ( i.LT.nlmnt) THEN
					DO j=i,nlmnt
						le(j)=le(j+1)
						lt(j)=lt(j+1)
						nd(j)=nd(j+1)
						np(j)=np(j+1)
						ln(1:nd(j),j)=ln(1:nd(j),j+1)
						lg(1:np(j),j)=lg(1:np(j),j+1)
					END DO
				END IF
			ELSE
				i= i + 1

			END IF

		END DO

!----------------------read nodal variables 
		DO i= 1, nnod
			i3= 3*i
			i2= i3 - 1
			i1= i2 - 1
			READ(in,ERR=991) ib(i)
			READ(in,ERR=991) x(i1), x(i2), x(i3),  dx(i1), dx(i2), dx(i3)
			READ(in,ERR=991) fc(i1), fc(i2), fc(i3)

		END DO

!---------------------read integration point variables
		DO i= 1, ngps
			READ(in,ERR=991) eps(i), dep(i)
			READ(in,ERR=991) ( sv(j,i), j= 1, nsv)
			READ(in,ERR=991) (dsv(j,i), j= 1, nsv)
			READ(in,ERR=991) (ep(j,i), j= 1, 6)
			READ(in,ERR=991) (de(j,i), j= 1, 6)
			READ(in,ERR=991) (st(j,i), j= 1, 6)
			READ(in,ERR=991) (ds(j,i), j= 1, 6)
			READ(in,ERR=991) (rt(j,i), j= 1, 3)
			READ(in,ERR=991) (dr(j,i), j= 1, 3)
		END DO

!----------------read (skip except e) material properties
		READ(in) nmatl
		DO i= 1, nmatl
			READ(in) nss(i)
			READ(in) twinlim
			READ(in) twinprob
		    	READ(in) e(1:9,i)
			
			DO j= 1, nss(i)
				READ(in)
			END DO

			DO j= 1, nss(i)
				READ(in) 
			END DO

			DO j=1,nss(i)
				READ(in) 
			END DO

		END DO

!---------------------------read (skip) surface properties
		READ(in) 

!---------------------------------read control data
		READ(in,ERR=991) ninc1, nxs, nfs
		DO i= 1, ninc1
			READ(in,ERR=991) tiem(i)
			IF ( nxs .GT. 0 ) THEN
				DO j= 1, nxs
					READ(in,ERR=991) ( x_spec(k,j,i), k= 1,15)
				END DO
			ENDIF
			IF ( nfs .GT. 0 ) THEN
				DO j= 1, nfs
					READ(in,ERR=991) ( f_spec(k,j,i), k= 1, 3)
				END DO
			ENDIF
		END DO

		l_file= .TRUE.						!  Set success flag to true


991	    CONTINUE							! -any failure by-passes flag 


	CLOSE(IN)


END SUBROUTINE FPREAD



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



SUBROUTINE 	INPFILE(in,l_file)
IMPLICIT 	none
LOGICAL(4) :: 		l_file
INTEGER ::			in,ierr
CHARACTER(100) ::	f_name,dir,dir1,dir2,dir3, dir4
CHARACTER(4) ::		ext

OPEN(in,FILE='dirname',STATUS='OLD', IOSTAT=ierr)
	IF ( ierr .EQ. 0 ) THEN
		READ(in,*) dir1,dir2,dir3,dir4
		dir=dir3
	ELSE
		WRITE(*,'("No dirname file found. Path is current directory.")')
		dir=''
	END IF
CLOSE (in)


WRITE(*,'("Name of FASOLT file (no extension): ")')
ext=".frs"
CALL check_file(f_name,ext,dir)
OPEN(in ,FILE=f_name ,STATUS='OLD', FORM= 'UNFORMATTED',	&
		&	ACCESS= 'SEQUENTIAL', IOSTAT=ierr)
	IF ( ierr .EQ. 0 ) THEN
		l_file = .TRUE.
		ELSE
		WRITE(*,'("Could not open results file! ")')
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

