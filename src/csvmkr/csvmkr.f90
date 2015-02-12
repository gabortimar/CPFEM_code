!***********************************************************************!
!																		!
!							     CVSMKR    								!
!							   Joao Fonseca                             !
!																		!
!-----------------------------------------------------------------------!
!           Post-pocessor: reads fst files and writes csv files 		!
!																		!
!																		!
!***********************************************************************!

PROGRAM CVSMKR
	
	IMPLICIT NONE
	INTEGER(4) type	
    
	CALL FPREAD()
	WRITE(*,'(/,10X,"Element (1,default) or Gauss point (2) values?")')
    	READ(*,*) type
	IF (type.EQ.2) THEN
		CALL GPSTATE()
	ELSE
    		CALL ELSTATE()
	END IF

	CALL POWER()

END PROGRAM



	
