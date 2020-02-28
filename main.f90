!MODULE sparse_system_solvers
!	IMPLICIT NONE
!	SAVE

!CONTAINS
!	FUNCTION solver_v1(A, b, m) result x
!		INTEGER, INTENT(in) 				:: m
!		REAL, DIMENSION(m,m), INTENT(in)	:: A
!		REAL, DIMENSION(m)					:: b
!		REAL, DIMENSION(m)					:: x
!	END FUNCTION solver_v1
!END MODULE

SUBROUTINE get_sparse_matrix(A, m)
	INTEGER, INTENT(in) 		:: m
	REAL, DIMENSION(m,m)		:: A
	REAL, DIMENSION(m)			:: rand_idx, rand_vect
	INTEGER, DIMENSION(m)		:: idx_x, idx_y
	INTEGER	i,j
	
	! Fill array elements with zeros
	A = 0
	
	! Initial diagonal elements with ones
	DO i=1,m
		A(i,i) = 1
	END DO

	! Take random m positions to initialize with random numbers
	CALL RANDOM_NUMBER(rand_idx)
	rand_idx = FLOOR((m+1)*rand_idx)
	idx_x = rand_idx

	CALL RANDOM_NUMBER(rand_idx)
	rand_idx = FLOOR((m+1)*rand_idx)
	idx_y = rand_idx

	CALL RANDOM_NUMBER(rand_vect)
	DO i=1,m
		A(idx_x(i),idx_y(i)) = rand_vect(i)
	END DO
END SUBROUTINE get_sparse_matrix


SUBROUTINE get_symmetric_sparse_matrix(A, m)
	INTEGER, INTENT(in) 		:: m
	REAL, DIMENSION(m,m)		:: A
	REAL, DIMENSION(m)			:: rand_idx, rand_vect
	INTEGER, DIMENSION(m)		:: idx_x, idx_y
	INTEGER	i,j
	
	! Fill array elements with zeros
	A = 0
	
	! Initial diagonal elements with ones
	DO i=1,m
		A(i,i) = 1
	END DO

	! Take random m positions to initialize with random numbers
	CALL RANDOM_NUMBER(rand_idx)
	rand_idx = FLOOR((m+1)*rand_idx)
	idx_x = rand_idx

	CALL RANDOM_NUMBER(rand_idx)
	rand_idx = FLOOR((m)*rand_idx)
	idx_y = rand_idx

	CALL RANDOM_NUMBER(rand_vect)
	DO i=1,m
		A(idx_x(i),idx_y(i)) = rand_vect(i)
		A(idx_y(i),idx_x(i)) = rand_vect(i)
	END DO
END SUBROUTINE get_symmetric_sparse_matrix





PROGRAM Main
	!USE sparse_system_solvers

	INTEGER, PARAMETER	::	m = 6		
	REAL		A(m,m)
	INTEGER i

	CALL get_sparse_matrix(A, m)
	print *, A
	print *, "------------"
	CALL get_symmetric_sparse_matrix(A, m)

	! COUT matrix to file
	DO i=1,m
		WRITE(12,*) A(i,:)
	END DO


END PROGRAM
