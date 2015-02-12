!***********************************************************************!
!																		!
!							  FASOLT 3									!
!							 Pete Bate									!
!-----------------------------------------------------------------------!
!						Three dimensional CPFEM							!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!

PROGRAM FASOLT

	USE FAS_COM			! FASOLT main variables declaration

	IMPLICIT NONE
	CHARACTER( LEN=4)	text
	CHARACTER( LEN=100) f_name
	INTEGER(4)			nsave, icurr, isave
	DIMENSION			isave(1:100)

!----------------------------------------------initialise and read data
	CALL FSINIT
	CALL RDDATA(f_name)
	CALL SETPTA

!----------------------------------------------get the output increments
	write(*,'(15X,''Number of increments is: '',I4)') ninc
	WRITE(*,'(/,10X,''Enter number of intermediate saves (max. 99):  '',\)')
	READ*,nsave
	IF( nsave .GT. 0) THEN
		WRITE(*,'(10X,''Enter intermediate increments to save at: '',\)')
		READ(*,*) ( isave(icurr),icurr=1,nsave)
	END IF
	
	! open file to write twinning info:
	open (unit=55,file="twins",action="write",status="replace")
!------------------------------------------perform incremental solutions
	DO inc= 1, ninc
		CALL UPDATE_NODE		!	update node etc. values
		CALL NRITER				!	iterate for solution
		CALL UPDATE_GPTV		!	update Gauss point values

		IF (tflag .ne. 0) THEN		! Check if there's any twinning.
			CALL DoReorient		! Do twin reorientation.
			tflag = 0		! reset overall twin flag to zero
		END IF

!-------------------------------!	save if required
		IF ( nsave .GT. 0 ) THEN
			DO icurr= 1, nsave
				IF ( inc .GE. isave(icurr) ) THEN
					WRITE( text, '(I4)') inc
					CALL WTDATA( TRIM(f_name)//'-'//TRIM(ADJUSTL(text))//'.frs')
					call WriteSlips( TRIM(f_name)//'-'//TRIM(ADJUSTL(text))//'.slp')
					isave(icurr)= ninc+1
				ENDIF
			END DO
		END IF

	END DO

!----------------------------------------write final result data
	close (55)
	WRITE( text, '(I4)') ninc
	CALL WTDATA( TRIM(f_name)//'-'//TRIM(ADJUSTL(text))//'.frs')
	call WriteSlips( TRIM(f_name)//'-'//TRIM(ADJUSTL(text))//'.slp')

END


!***********************************************************************!
!																		!
!         SUBROUTINE DoReorient						!
!-----------------------------------------------------------------------!
!	Reorient twinned IPs according to angles specified in fmt file				!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE DoReorient

	USE FAS_COM			! FASOLT main variables declaration
	USE FAS_CPL			! FASOLT material property block

	IMPLICIT NONE

	INTEGER(4) 	i, j
	REAL(4) 	e1, e2, e3, e1n, e2n, e3n, c1, c2, c3, ss1, s2, s3, ang, c, s, u1, u2, u3
			! ss1 to avoid conflict with module

	REAL(4), DIMENSION(1:3)::	nvec, dvec
	REAL(4), DIMENSION(1:3, 1:3)::	RM, ERM, FRM


	 do j = 1, ngps				! loop over IPs
	  if (twinsys(j) .ne. 0) then	! reorient this IP!
		! get index of active twinning system:
		i = twinsys(j)
		
		! get current Euler angles and twin reorientation angle:
		e1 = sv(1,j)
		e2 = sv(2,j)
		e3 = sv(3,j)
		ang = ssang(i,1) * 3.14159 / 180.	! convert to radians

		! calculate rotation matrix of IP from current Euler angles:
		c1 = COS(e1)
		c2 = COS(e2)
		c3 = COS(e3)
		ss1 = SIN(e1)
		s2 = SIN(e2)
		s3 = SIN(e3)

		RM(1,1) = c1*c3 - c2*ss1*s3
		RM(1,2) = -c1*s3 - c2*c3*ss1
		RM(1,3) = ss1*s2
		RM(2,1) = c3*ss1 + c1*c2*s3
		RM(2,2) = c1*c2*c3 - ss1*s3
		RM(2,3) = -c1*s2
		RM(3,1) = s2*s3
		RM(3,2) = c3*s2
		RM(3,3) = c2

		! calculate current normal vector and direction vector:
		nvec = MATMUL(RM, an(1:3,i,1))
		dvec = MATMUL(RM, ab(1:3,i,1))

		! calculate rotation axis for twin reorientation:
		u1 = nvec(2) * dvec(3) - nvec(3) * dvec(2)
		u2 = nvec(3) * dvec(1) - nvec(1) * dvec(3)
		u3 = nvec(1) * dvec(2) - nvec(2) * dvec(1)

		! calculate extra rotation matrix for twin reorientation:
		c = COS(ang)
		s = SIN(ang)		

		ERM(1,1) = c + u1*u1 * (1. - c)
		ERM(1,2) = u1 * u2 * (1. - c) - u3 * s
		ERM(1,3) = u1 * u3 * (1. - c) + u2 * s
		ERM(2,1) = u1 * u2 * (1. - c) + u3 * s
		ERM(2,2) = c + u2*u2 * (1. - c)
		ERM(2,3) = u2 * u3 * (1. - c) - u1 * s
		ERM(3,1) = u1 * u3 * (1. - c) - u2 * s
		ERM(3,2) = u2 * u3 * (1. - c) + u1 * s
		ERM(3,3) = c + u3*u3 * (1. - c)

		! multiply matrices to get updated rotation matrix for IP:
		FRM = MATMUL(ERM, RM)

		! calculate new Euler angles from updated rotation matrix:
		e1n = ATAN2( FRM(1,3) , -FRM(2,3) )
		! check bounds for FRM(3,3) - if out of bounds, normalize:
		if (abs(FRM(3,3)) .gt. 1.) FRM(3,3) = FRM(3,3) / abs(FRM(3,3))
		e2n = ACOS( FRM(3,3) )
		e3n = ATAN2( FRM(3,1) , FRM(3,2) )

		! update sv variable using new Euler angle values:
		sv(1,j) = e1n
		sv(2,j) = e2n
		sv(3,j) = e3n
		! IP reorientation done!

		! write twinning information to file:
		! (incr. number, IP index, system index, Eulers before, Eulers after)
		write(55,'(I8, I8, I8, 9999ES15.5)') inc, j, i, e1, e2, e3, e1n, e2n, e3n

	  	twinsys(j) = 0	! reset twinning system to 0!
		! increment number of times this IP has twinned:
		timestwinned(j) = timestwinned(j) + 1
	  end if

	 end do






END SUBROUTINE DoReorient

!***********************************************************************!
!																		!
!         SUBROUTINE FSINIT												!
!-----------------------------------------------------------------------!
!			Initialises FASOLT 3 and gets control data					!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE FSINIT

	USE FAS_COM			! FASOLT main variables declaration

	IMPLICIT NONE

	INTEGER(4)	iomon, in, out, ierr

    COMMON/ iochan/ iomon, in, out

	CHARACTER(100) :: dir,dir1,dir2,dir3,dir4,fname
!----------------------------------------integer channel numbers
	iomon= 51
	in=    52
	out=   53

	OPEN(in,FILE='dirname',STATUS='OLD', IOSTAT=ierr)
		IF ( ierr .EQ. 0 ) THEN
			READ(in,*) dir1,dir2,dir3,dir4,dir
		ELSE
			WRITE(*,'("No dirname file found. Path is current directory.")')
			dir='./'
		END IF
	CLOSE(in)

	fname=TRIM(dir)//'fafcdf.dat'
	
	OPEN( UNIT=in, FILE=fname)
       READ( in, *) toler, maxit, alpha, afactor
       READ( in, *) efix, smalp, smalf
       READ( in, *) smald, smale, smalr
	CLOSE( in)

END SUBROUTINE FSINIT

!***********************************************************************!
!																		!
!         SUBROUTINE UPDATE-NODE										!
!-----------------------------------------------------------------------!
!		Updates relevant variables at START of incremental solution.	!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE UPDATE_NODE

	USE FAS_COM			! FASOLT main variables declaration

	IMPLICIT NONE
	INTEGER(4)	ndof

	ndof= 3*nnod
	x(1:ndof)=x(1:ndof) + dx(1:ndof)		!	update node positions

!---------------------------------------Update die positions


END SUBROUTINE UPDATE_NODE

!***********************************************************************!
!																		!
!         SUBROUTINE UPDATE-GPTV										!
!-----------------------------------------------------------------------!
!		Updates relevant variables at END of incremental solution.		!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE UPDATE_GPTV

	USE FAS_COM			! FASOLT main variables declaration
	
	IMPLICIT NONE

	INTEGER(4)	i, j

!-------------------------------------update all Gauss point values
	DO  i= 1, ngps

		eps( i)= eps( i) + dep( i)
		DO j=1,3
			rt( j, i)= rt( j, i) + dr( j, i)	! update rotation vector
		END DO
			
		DO j=1,6
			ep( j, i)= ep( j, i) + de( j, i)
			st( j, i)= st( j, i) + ds( j, i)
		END DO

        	DO j=1,nsv
			sv( j, i)= sv( j, i) + dsv( j, i)	! update Euler angles (+ slip resistances)
		END DO

	END DO

!--------------------------------dump selected Gauss point values


END SUBROUTINE UPDATE_GPTV

!***********************************************************************!
!																		!
!         SUBROUTINE FAFERR(IER)										!
!-----------------------------------------------------------------------!
!        -Error reporting routine with possible dump and stop			!
!-----------------------------------------------------------------------!
!         Arguments       ier   INTEGER, Error code:					!
!                                               ier>500 fatal			!
!                                               ier<500 warning			!
!                                               ier>900 fatal, dump		!
!-----------------------------------------------------------------------!
!																		!
!***********************************************************************!
SUBROUTINE FAFERR(ier)

	IMPLICIT NONE
	
	INTEGER(4)		iomon, in, out, ier

    COMMON/ iochan/ iomon, in, out

	SELECT CASE (ier)

		CASE (:500)
			WRITE(iomon, '(/,10X,'' *** FAFNER:  WARNING, '',I3,/)')  ier
			WRITE(*    , '(/,10X,'' *** FAFNER:  WARNING, '',I3,/)')  ier

		CASE (501:)
			WRITE(iomon, '(/,10X,'' *** FAFNER:  FATAL ERROR, '',I3,/)') ier
			WRITE(*    , '(/,10X,'' *** FAFNER:  FATAL ERROR, '',I3,/)') ier

			IF ( ier .GT. 900 ) CALL WTDATA('faferr.frs')
			WRITE(*,'(10x,''any key to finish..'',\)')
			READ(*,*)
							STOP
     END SELECT

END SUBROUTINE FAFERR

