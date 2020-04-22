MODULE timings
    IMPLICIT NONE
    SAVE
    PRIVATE
    PUBLIC tic,toc,startClock,stopClock
    INTEGER, PARAMETER 		:: dp = selected_real_kind(15,307)
    INTEGER, PARAMETER 		:: dp_int = selected_int_kind(15)
	REAL(KIND=dp) 			:: startT, finishT
	REAL(KIND=dp)		 	:: count_rate
	INTEGER(KIND=dp_int) 	:: countS, countF

CONTAINS
    SUBROUTINE tic(startTime)
		REAL(KIND=dp), optional :: startTime
		IF (present(startTime)) THEN
			CALL cpu_time(startTime)
		ELSE
			CALL cpu_time(startT)
		ENDIF
    END SUBROUTINE tic


    SUBROUTINE toc(elapsedTime, startTime)
		REAL(KIND=dp), optional :: elapsedTime, startTime
		CALL cpu_time(finishT)
		IF (present(elapsedTime) .AND. present(startTime)) THEN
			elapsedTime = finishT - startTime
			PRINT *, "ELAPSED CPU TIME SINCE tic():", elapsedTime
		ELSEIF (present(elapsedTime)) THEN
			elapsedTime = finishT - startT
		ELSE
			PRINT *, "ELAPSED CPU TIME SINCE tic():", finishT - startT
		ENDIF
    END SUBROUTINE toc

    SUBROUTINE startClock(startTime)
		INTEGER(KIND=dp_int), optional 	:: startTime

		IF (present(startTime)) THEN
			CALL SYSTEM_CLOCK(startTime, count_rate)
		ELSE
			CALL SYSTEM_CLOCK(countS, count_rate)
		ENDIF
    END SUBROUTINE startClock
    
    SUBROUTINE stopClock(elapsedTime, startTime)
		REAL(KIND=dp), optional 		:: elapsedTime
		INTEGER(KIND=dp_int), optional 	:: startTime
		CALL SYSTEM_CLOCK(countF)

		IF (present(elapsedTime) .AND. present(startTime)) THEN
			elapsedTime = (countF - startTime)/count_rate
			PRINT *, "ELAPSED WALL CLOCK TIME:", elapsedTime
		ELSEIF (present(elapsedTime)) THEN
			elapsedTime = (countF - countS)/count_rate
		ELSE
			PRINT *, "ELAPSED WALL CLOCK TIME:", (countF - countS)/count_rate
		ENDIF
    END SUBROUTINE stopClock
END
