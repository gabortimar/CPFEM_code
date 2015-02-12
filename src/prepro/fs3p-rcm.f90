!*******************************************************************************!
!																				!
!		SUBROUTINE RENUMBER														!
!        reverse Cuthill-McGee numbering optimisation							!
!																				!
!*******************************************************************************!
SUBROUTINE RENUMBER(iprof)
      
	IMPLICIT NONE
	LOGICAL(4)	again, swap
	INTEGER(4)	maxlev, lsv, lcs, isv, i, j, k, lct, ibt, nm, minc,		&
										&	level, maxl, ii, iprof, lnot
	REAL(4)		xt

    DIMENSION lsv(1:20000),lcs(1:20000),isv(1:20000)

!----------------------------------FAFNER 3-D common blocks and declarations
	INTEGER(4)	nlmnt, nnod, ngps, ninc, ninc0, ninc1, le, lt, nd, np, ln, lg, ib
	REAL(4)		x, dx, fc

	COMMON/ numbrs/ nlmnt, nnod, ngps, ninc, ninc0, ninc1
    COMMON/ elnums/ le(1000), lt(1000), nd(1000), np(1000), ln(20,1000), lg(8,1000)
    COMMON/ glbal1/ ib(5000), x(15000), dx(15000), fc(15000)


    DATA maxlev/20000/		! maximum number of levels


	IF ( nlmnt .LE. 1 )		RETURN		!	don't even try

!----------------------initialise 'number connected-to-node' list
	DO i= 1, nnod
		IF( ib(i) .EQ. 0 )THEN
			lcs(i)= 0
		ELSE
			lcs(i)= 1
		ENDIF
	END DO

    DO i= 1, nlmnt
		DO j= 1, nd(i)
			k= ln(j,i)
			lcs(k)= lcs(k)+1
		END DO
	END DO

!------------------------get start and generate trial structure
    nm= 0
    minc= 10
	DO i= 1, nnod
		IF ( ( lcs(i) .LT. minc ) .AND. ( lcs(i) .GT. 2 ) ) THEN
			minc= lcs(i)
			nm= i
		ENDIF
	END DO

    CALL ROLSTRUC(.FALSE.,maxlev,nm,lsv,lcs,isv,level)
    maxl=level

!-------------------------------search for maximum graph diameter
	again=.TRUE.
	do while (again)

		again=.FALSE.
		DO  i= isv(level), isv(level-1)+1, -1
			ii=lsv(i)
			CALL ROLSTRUC(.FALSE.,maxlev,ii,lsv,lcs,isv,level)

			IF ( level .GT. maxl ) THEN  !better one found
				again=.TRUE.
				maxl=level
				nm= ii
			ENDIF

		END DO

	END DO

    CALL ROLSTRUC(.TRUE.,maxlev,nm,lsv,lcs,isv,iprof) !	re-generate best structure

!------------------------------- reverse-number re-order nodal sequence
    DO i= 1, nlmnt
		DO j= 1, nd(i)
			k= ln(j,i)
			ln(j,i)= lsv(k)
		END DO
	END DO

	again=.true.
    lnot= nnod + 1
	DO WHILE ( again)
		lnot=lnot-1

		swap=.FALSE.
		DO i= 2, lnot
			IF( lsv(i) .LT. lsv(i-1) )THEN
				swap=.TRUE.

				lct= lsv(i)
				lsv(i)= lsv(i-1)
				lsv(i-1)= lct

				DO j= -2,0
					xt= x(3*i +j)
					x(3*i +j)= x(3*(i-1) +j)
					x(3*(i-1) +j)=xt
				END DO	

				ibt= ib(i)
				ib(i)= ib(i-1)
				ib(i-1)= ibt

			ENDIF
		END DO

		AGAIN= ( swap .AND. (lnot .GT. 2) )
	END DO

END SUBROUTINE

!*******************************************************************************!
!																				!
!		SUBROUTINE ROLSTRUC														!
!              Produces a rooted, ordered level structure						!
!																				!
!*******************************************************************************!
SUBROUTINE ROLSTRUC(renum,nw,iroot,lsv,lcs,isv,iprof)

	IMPLICIT NONE
    LOGICAL(4)	renum, conx
	INTEGER(4)	nw, iroot, lsv, lcs, isv, iprof, i, i1, il, in, j, jm, jmin, &
							&	jnn, jv, jvm, k, l, num, nums, level, lst, lfn

    DIMENSION	lsv(1:*),lcs(1:*),isv(1:*)

!----------------------------------FAFNER 3-D common blocks and declarations
	INTEGER(4)	nlmnt, nnod, ngps, le, lt, nd, np, ln, lg, ninc, ninc0, ninc1

	COMMON/ numbrs/ nlmnt, nnod, ngps, ninc, ninc0, ninc1
    COMMON/ elnums/ le(1000), lt(1000), nd(1000), np(1000), ln(20,1000), lg(8,1000)


!-----------------------------------------sort out the structure by magic
	DO i=1,nnod
		lcs(i)=-lcs(i)
	END DO

    lsv(1)= iroot
    lcs(iroot)= -lcs(iroot)
    isv(1)= 1
    num= 1
    level= 1

	DO WHILE ( num .LT. nnod )

		IF( level .LT. 2 ) THEN
			lst= 1
		ELSE
			lst= isv(level-1)+1
		ENDIF

		lfn= isv(level)
		level= level +1

		DO i= lst, lfn
			in= lsv(i)
			nums= num +1

			DO il= 1, nlmnt
				conx= .FALSE.
				jnn= nd(il)
				DO j= 1, jnn
					IF( ln(j,il) .EQ. in )  conx=.TRUE.
				END DO

				IF ( conx ) THEN
					DO j=1,jnn
						l= ln(j,il)

						IF ( lcs(l) .LT. 0 ) THEN
							num= num +1
							lcs(l)= -lcs(l)
							k= num

							IF ( num .GT. nums ) THEN

								DO i1= num-1, nums, -1
									IF ( lcs(l) .LT. lcs(lsv(i1)) ) k= i1
								END DO

								IF ( k .LT. num ) THEN
									DO i1= num, k+1, -1
										lsv(i1)= lsv(i1-1)
									END DO
								ENDIF

							ENDIF
							lsv(k)=l

						ENDIF

					END DO
				ENDIF
			END DO
		END DO
		isv(level)= num

    END DO


	IF( renum ) THEN	! form re-ordering sequence and calculate profile

		DO i= 1, nnod
			j= nnod -i +1
			isv(lsv(i))= j
		END DO

		DO i= 1, nnod
			lsv(i)= isv(i)
			isv(i)= 1
		END DO

		DO i=1,nlmnt
			jm= nd(i)
			jmin=lsv(ln(1,i))

            DO j= 2, jm
				jv= lsv(ln(j,i))
				IF( jv .LT. jmin )  jmin= jv
			END DO

			DO j=1,jm
				jv= lsv(ln(j,i))
				jvm= jv- jmin +1
				IF( jvm .GT. isv(jv) )  isv(jv)= jvm
			END DO
		END DO

		iprof= 0
		DO i= 1, nnod
			iprof=iprof +isv(i)
		END DO

	ELSE				!just return number of levels as profile size
		iprof= level

	ENDIF

END SUBROUTINE


