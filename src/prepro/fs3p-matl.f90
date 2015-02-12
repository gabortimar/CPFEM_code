!*******************************************************************************!
!																				!
!	SUBROUTINE SETMATC															!
!		Sets the materials parameters											!
!																				!
!																				!
!*******************************************************************************!
SUBROUTINE SETMATC

	USE FAS_COM

	IMPLICIT NONE
	CHARACTER(LEN=1) cdum
	LOGICAL(4)	l_file
	INTEGER(4)	i, j, k, in, i_plas
	REAL(4)		inislr, sumf, r, dummy, f_elas

	REAL(4), DIMENSION(1:maxss)::	ini_ss
	REAL(4), DIMENSION(1:12)::		frac

!---------------------------------------------------no. of property sets to use
	WRITE(*,'(20X,''Material Property Data'',/)')
	WRITE(*,'(20X,''Enter number of property sets (max.12)'',\)')
	READ(*,*) nmatl


!------------------------------------------get the number of element codes set
	j=0
	DO i= 1, nlmnt
		IF (le(i) .GT. j) j=le(i)
	END DO

	IF ( j .GT. nmatl ) THEN
		WRITE(*,'(/,15X,''There are '',I3,'' element codes! '')') j
		WRITE(*,'(/,15X,''Assigning at random: enter fraction for: '')')
		DO i= 1, nmatl
			WRITE(*,'(18X,''Material '',I3,'': '',\)') i
			READ(*,*) frac(i)
		END DO

!------------------------------normalise and make cummulative
		sumf= SUM(frac(1:nmatl)) ; frac(1:nmatl)= frac(1:nmatl)/sumf
		IF ( nmatl.GT.1) THEN
			DO i=2, nmatl
				frac(i)= frac(i)+frac(i-1)
			END DO
		END IF

		DO i= 1, nlmnt
			CALL RANDOM_NUMBER(r)
			k=0
			DO j= 1, nmatl
				IF ( r .GT.frac(j) ) k=j
			END DO
			le(i)= k+1	
		END DO
	END IF

!---------------------------------------------get the material property data
	nmpar=0
	ampar=0.

	DO mat_code = 1, nmatl
		WRITE(*,'(25X,''<enter> to input property set '',I3,\)'),mat_code
		READ(*,*) 
		
		IN=51
		CALL MATFILE(in,l_file)

		READ(in,*) i_plas

		IF (mat_code.EQ.1) THEN
			nmpar= i_plas
		ELSE
			IF (i_plas.NE.nmpar) THEN
				WRITE(*,'(10X,''Inconsistent plastic property list!'')')
				PAUSE
				STOP
			END IF
		END IF

		READ(in,*) nss(mat_code)

!------------------------------------read elastic constants
		READ(in,*) (e(i, mat_code), i= 1, 9)

!------------------------------------read ss angles (also indicator of twin)
		READ(in,*) (ssang(i, mat_code), i= 1, nss(mat_code))

!-----------------------------------------read slip systems
		DO i= 1, nss(mat_code)
			READ(in,*) (ab(j,i,mat_code),j=1,3), &
					&  (an(j,i,mat_code),j=1,3)
		END DO
!---------------------------------------read ss properties
		DO i= 1, nss(mat_code)
			READ(in,*) (ampar(j,i,mat_code),j=1,nmpar)
		END DO

!--------------------------------read latent hardening matrix
		DO i=1,nss(mat_code)
			READ(in,*) (smatrx(j,i, mat_code),j=1,nss(mat_code))
		END DO
        CLOSE(IN)
!---------------------------------------------get the initial slip resistances
		WRITE(*,'(5X,''Make initial slip resistances equal? (Y/N <ret>:'',\)')
		READ(*,*) cdum

		IF ( cdum.EQ.'y' .OR. cdum.EQ.'Y') THEN
			WRITE(*,'(5X,''Enter the initial slip resistance: '',\)')
			READ(*,*) inislr
			ini_ss(1:nss(mat_code))= inislr
		ELSE
			WRITE(*,'(5X,''Enter the '',I3,	&
				&	'' initial slip resistances:'',\)'),nss(mat_code)
			READ(*,*) ini_ss(1:nss(mat_code))
		END IF

		DO i=1,nlmnt
			IF ( le(i) .EQ. mat_code ) THEN
				DO j=1,np(i)
					sv(4:3+nss(mat_code),lg(j,i)) = ini_ss(1:nss(mat_code))
				END DO
			END IF
		END DO

	END DO

!--------------------------------------------------option for random elastic fraction
	WRITE(*,'(/,2X,''Set a fraction to elastic only? (y/n): '',\)') 
	READ(*,*) cdum
	IF ( cdum.EQ.'y' .OR. cdum.EQ.'Y' ) THEN
			WRITE(*,'(7X,''Enter elastic fraction : '',\)') 
			READ(*,*) f_elas

			j=0
			DO i= 1, nlmnt

				CALL RANDOM_NUMBER(dummy)
				IF ( dummy .LE. f_elas ) THEN
					le(i) =-le(i) ; j=j+1
				END IF

			END DO

			WRITE(*,'(7X, G12.5,'' set elastic '')') REAL(j, 4)/ REAL(nlmnt, 4)

	END IF


END SUBROUTINE

