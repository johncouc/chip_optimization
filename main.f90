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




			!REAL(selected_real_kind(6,37)), DIMENSION(:,:), INTENT(in)				:: A
			!REAL(selected_real_kind(6,37)), DIMENSION(SIZE(A,DIM=1), SIZE(A,DIM=2))	:: AA
			!REAL(selected_real_kind(6,37)), DIMENSION(:),  INTENT(out)				:: lambdar, lambdai
			!REAL(selected_real_kind(6,37)), DIMENSION(:,:), ALLOCATABLE				:: VLR
		  	!REAL(selected_real_kind(6,37)), DIMENSION(:), ALLOCATABLE				:: work
			!INTEGER info, l_work
			!AA = A
			!ALLOCATE(work(1))
			!CALL SGEEV('N','N',SIZE(A, DIM=2), AA, SIZE(A, DIM=2), lambdar, lambdai, VLR, 1, VLR, 1, work, -1, info)
			!l_work = work(1)
			!DEALLOCATE(work)
			!ALLOCATE(work(l_work))
			!CALL SGEEV('N','N',SIZE(A, DIM=2), AA, SIZE(A, DIM=2), lambdar, lambdai, VLR, 1, VLR, 1, work, l_work, info)
			!DEALLOCATE(work)			



	END SUBROUTINE solver_v1
END MODULE sparse_system_solvers

SUBROUTINE get_sparse_matrix(A, m, rand_idx, rand_vect_real)
	INTEGER, INTENT(in) 				:: m
	INTEGER, DIMENSION(m,m)				:: A
	REAL, DIMENSION(m), INTENT(IN)		:: rand_vect_real
	REAL, DIMENSION(2*m), INTENT(IN)	:: rand_idx
	REAL, DIMENSION(2*m)				:: rand_idx_COPY
	INTEGER, DIMENSION(m)				:: idx_x, idx_y, rand_vect_int
	INTEGER	i,j
	
	! Fill array elements with zeros
	A = 0
	
	! Initial diagonal elements with ones
	DO i=1,m
		A(i,i) = 1
	END DO

	rand_idx_COPY = rand_idx

	! SET RANDOM INDICES TO INTEGERS
	rand_idx_COPY = FLOOR((m)*rand_idx_COPY)
	DO i=1,m
		idx_x(i) = rand_idx_COPY(i)
		idx_y(i) = rand_idx_COPY(i+m)
	END DO

	! SET RANDOM VALUES TO INTEGERS
	rand_vect_int = FLOOR((m)*rand_vect_real)
	DO i=1,m
		A(idx_x(i),idx_y(i)) = rand_vect_int(i)
	END DO
END SUBROUTINE get_sparse_matrix


SUBROUTINE get_symmetric_sparse_matrix(A, m, rand_idx, rand_vect_real)
	INTEGER, INTENT(in) 				:: m
	INTEGER, DIMENSION(m,m)				:: A
	REAL, DIMENSION(m)					:: rand_vect_real
	REAL, DIMENSION(2*m), INTENT(IN)	:: rand_idx
	REAL, DIMENSION(2*m)				:: rand_idx_COPY
	INTEGER, DIMENSION(m)				:: idx_x, idx_y, rand_vect_int
	INTEGER	i,j
	
	! Fill array elements with zeros
	A = 0
	
	! Initial diagonal elements with ones
	DO i=1,m
		A(i,i) = 1
	END DO

	rand_idx_COPY = rand_idx

	! SET RANDOM INDICES TO INTEGERS
	rand_idx_COPY = FLOOR((m)*rand_idx_COPY)
	DO i=1,m
		idx_x(i) = rand_idx_COPY(i)
		idx_y(i) = rand_idx_COPY(i+m)
	END DO
	
	! SET RANDOM VALUES TO INTEGERS
	rand_vect_int = FLOOR((m)*rand_vect_real)
	DO i=1,m
		A(idx_x(i),idx_y(i)) = rand_vect_int(i)
		A(idx_y(i),idx_x(i)) = rand_vect_int(i)
	END DO

END SUBROUTINE get_symmetric_sparse_matrix

SUBROUTINE get_random_vector(b, m)
	INTEGER, INTENT(in) 	:: m
	REAL, DIMENSION(m)		:: b

	CALL RANDOM_SEED()
	CALL RANDOM_NUMBER(b)
END SUBROUTINE get_random_vector

PROGRAM Main
	!USE sparse_system_solvers

	INTEGER, PARAMETER		:: m = 6		
	INTEGER, DIMENSION(m,m)	:: A
	REAL, DIMENSION(m)		:: b, rand_vect_real
	REAL, DIMENSION(2*m)	:: rand_idx
	INTEGER i,n

	! RANDOM MATRIX
	CALL RANDOM_SEED()
	CALL RANDOM_NUMBER(rand_idx)

	CALL RANDOM_SEED()
	CALL RANDOM_NUMBER(rand_vect_real)

	CALL get_sparse_matrix(A, m, rand_idx, rand_vect_real)
	CALL get_symmetric_sparse_matrix(A, m, rand_idx, rand_vect_real)

	! COUT matrix to file
	DO i=1,m
		WRITE(12,*) A(i,:)
	END DO
END PROGRAM
