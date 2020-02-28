MODULE sparse_system_solvers
	IMPLICIT NONE
	SAVE

CONTAINS
	SUBROUTINE solver_v1(A, b, m, x)
		INTEGER, INTENT(in) 				:: m
		REAL, DIMENSION(m,m), INTENT(in)	:: A
		REAL, DIMENSION(m)					:: b
		REAL, DIMENSION(m), INTENT(out)		:: x
		PRINT *, "AHA"




			REAL(selected_real_kind(6,37)), DIMENSION(:,:), INTENT(in)				:: A
			REAL(selected_real_kind(6,37)), DIMENSION(SIZE(A,DIM=1), SIZE(A,DIM=2))	:: AA
			REAL(selected_real_kind(6,37)), DIMENSION(:),  INTENT(out)				:: lambdar, lambdai
			REAL(selected_real_kind(6,37)), DIMENSION(:,:), ALLOCATABLE				:: VLR
		  	REAL(selected_real_kind(6,37)), DIMENSION(:), ALLOCATABLE				:: work
			INTEGER info, l_work
			AA = A
			ALLOCATE(work(1))
			CALL SGEEV('N','N',SIZE(A, DIM=2), AA, SIZE(A, DIM=2), lambdar, lambdai, VLR, 1, VLR, 1, work, -1, info)
			l_work = work(1)
			DEALLOCATE(work)
			ALLOCATE(work(l_work))
			CALL SGEEV('N','N',SIZE(A, DIM=2), AA, SIZE(A, DIM=2), lambdar, lambdai, VLR, 1, VLR, 1, work, l_work, info)
			DEALLOCATE(work)			



	END SUBROUTINE solver_v1
END MODULE sparse_system_solvers

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

SUBROUTINE get_random_vector(b, m)
	INTEGER, INTENT(in) 	:: m
	REAL, DIMENSION(m)		:: b

	CALL RANDOM_NUMBER(b)
END SUBROUTINE get_random_vector





PROGRAM Main
	!USE sparse_system_solvers

	INTEGER, PARAMETER		:: m = 6		
	REAL, DIMENSION(m,m)	:: A
	REAL, DIMENSION(m)		:: b
	INTEGER i

	CALL get_sparse_matrix(A, m)
	print *, A
	print *, "------------"
	CALL get_symmetric_sparse_matrix(A, m)

	
	CALL get_random_vector(b, m)

	! COUT matrix to file
	DO i=1,m
		WRITE(12,*) A(i,:)
	END DO


END PROGRAM
