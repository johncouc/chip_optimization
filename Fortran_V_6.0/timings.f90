module timings
    implicit none
    save
    private
    public tic,toc,startClock,stopClock
    integer, parameter 		:: dp = selected_real_kind(15,307)
    integer, parameter 		:: dp_int = selected_int_kind(15)
	real(kind=dp) 			:: startT, finishT
	real(kind=dp)		 	:: count_rate
	integer(kind=dp_int) 	:: countS, countF

contains
    subroutine tic(startTime)
		real(kind=dp), optional :: startTime
		if (present(startTime)) THEN
			call cpu_time(startTime)
		else
			call cpu_time(startT)
		endif
    end subroutine tic


    subroutine toc(elapsedTime, startTime)
		real(kind=dp), optional :: elapsedTime, startTime
		call cpu_time(finishT)
		if (present(elapsedTime) .and. present(startTime)) then
			elapsedTime = finishT - startTime
			print *, "ELAPSED CPU TIME SINCE tic():", elapsedTime
		elseif (present(elapsedTime)) then
			elapsedTime = finishT - startT
		else
			print *, "ELAPSED CPU TIME SINCE tic():", finishT - startT
		endif
    end subroutine toc

    subroutine startClock(startTime)
		integer(kind=dp_int), optional 	:: startTime

		if (present(startTime)) THEN
			CALL SYSTEM_CLOCK(startTime, count_rate)
		else
			call SYSTEM_CLOCK(countS, count_rate)
		endif
    end subroutine startClock
    
    subroutine stopClock(elapsedTime, startTime)
		real(kind=dp), optional 		:: elapsedTime
		integer(kind=dp_int), optional 	:: startTime
		call SYSTEM_CLOCK(countF)

		if (present(elapsedTime) .and. present(startTime)) then
			elapsedTime = (countF - startTime)/count_rate
			print *, "ELAPSED WALL CLOCK TIME:", elapsedTime
		elseif (present(elapsedTime)) then
			elapsedTime = (countF - countS)/count_rate
		else
			print *, "ELAPSED WALL CLOCK TIME:", (countF - countS)/count_rate
		endif
    end subroutine stopClock
end
