! COMMANDS TO RUN:
! gfortran -c main.f90 eigv_module.f90
! gfortran -c main.f90 eigv_module.f90
! ./main
! The matrices and results will be saved in output files 'fort.*'
! where * respresents 11,12 for the input matrix and vector
! 9 and 10 - output solutions for two methods

! Module is based on: https://www.netlib.org/lapack/lug/node38.html
MODULE sparse_system_solvers
	IMPLICIT NONE
	SAVE

CONTAINS

	! Solves general band matrices including tridiagonal matrices:
	! LU factorization with partial pivoting
	SUBROUTINE solver_v1(A, b, m, x)
		INTEGER								  	:: m
		REAL, DIMENSION(m,m), INTENT(in)		:: A
		REAL, DIMENSION(m,m)					:: AA
		REAL, DIMENSION(m), INTENT(in)			:: b
		REAL, DIMENSION(m), INTENT(out)			:: x
		INTEGER, DIMENSION(m)					:: IPIV
		INTEGER info, l_work
		
		! Make matrix and vector copies
		AA = A
		x = b

		! Compute LU factorization 
		CALL SGETRF(SIZE(AA,1), SIZE(AA,2), AA, m, IPIV, info)

		! Solve system
		CALL SGETRS('N',SIZE(AA, DIM=2),1, AA, SIZE(AA, DIM=1), IPIV, x, SIZE(x,1), info)
	END SUBROUTINE solver_v1


	! Solves symmetric and positive definite matrices including band matrices:
	! Cholesky factorization
	SUBROUTINE solver_v2(A, b, m, x)
		INTEGER								  	:: m
		REAL, DIMENSION(m,m), INTENT(in)		:: A
		REAL, DIMENSION(m,m)					:: AA, A_u
		REAL, DIMENSION(m), INTENT(in)			:: b
		REAL, DIMENSION(m), INTENT(out)			:: x
		INTEGER, DIMENSION(m)					:: IPIV
		INTEGER info, l_work
		
		! Make matrix and vector copies
		AA = A
		x = b

		! Solves a system of linear equations A*X = B where 
		! A is an N-by-N symmetric positive definite matrix
		! Performs Cholesky factorization
		CALL SPOSV('L', SIZE(AA,1), 1, AA, SIZE(AA,1), x, SIZE(x,1), info)
	END SUBROUTINE solver_v2

	! BASED ON THIS COMMENT: https://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=5048
	! TRY BAND SOLVERS IN LLAPACK
	! Computes the solution to a real system of linear equations A * X = B,
 	! where A is an N-by-N symmetric positive definite band matrix and X
 	! and B are N-by-NRHS matrices.
	SUBROUTINE solver_v3(A, b, m, x)
		INTEGER								  	:: m, nsuper, info, i, j
		REAL, DIMENSION(m,m), INTENT(in)		:: A
		REAL, DIMENSION(m,m)					:: A_u
		REAL, DIMENSION(CEILING(m/2.0),m)			:: A_banded
		REAL, DIMENSION(m)						:: b, x

		x = b
		nsuper = 1
		A_u = 0.0
		A_banded = 0.0

		! GET THE UPPER DIAGONAL PART
        CALL SLACPY('U', m, m, A, m, A_u, m)

		DO j=1,m
			DO i=1,m
				IF (MAX(1,j-nsuper)<=i .AND. i<=j) THEN
					A_banded(nsuper+1+i-j,j) = A_u(i,j)	
				END IF
			END DO
		END DO
		
		CALL SPBSV('U', m, nsuper, 1, A_banded, size(A_banded, 1), x, m, info)		
		PRINT *, "Is system solved? (0-yes):  ", info
	END SUBROUTINE solver_v3
END MODULE sparse_system_solvers

! DOES NOT WORK: IDEA - set work variable properly - allocate, deallocate
SUBROUTINE get_sparse_matrix_llapack(A, m)
	INTEGER, INTENT(in) 				:: m
	REAL, DIMENSION(m,m)				:: A
	INTEGER, DIMENSION(4)				:: iseed
	REAL, DIMENSION(m)					:: diag, dl, dr, ipivot
	INTEGER, DIMENSION(m)				:: iwork
	REAL sparsity, anorm
	INTEGER info


	! Define the seed
	iseed = 4

	! Define diagonal, sparsity percentage
	diag = 1
	sparsity = 0.8
	anorm = 5.0 ! MAX ENTRY OF OUTPUT MATRIX

	!CALL SLATMR(m,m,'N',iseed,'S',diag,0,1.0,1.0,'F','N',dl,0,0.0,dr,0,0.0,'N',ipivot,sparsity,0,0,anorm,'N',A,m,iwork,info)
END SUBROUTINE get_sparse_matrix_llapack


SUBROUTINE get_sparse_matrix(A, m, rand_idx, rand_vect_real)
	INTEGER, INTENT(in) 				:: m
	REAL, DIMENSION(m,m)				:: A
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
	REAL, DIMENSION(m,m)				:: A
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
	INTEGER				 	:: m
	REAL, DIMENSION(m)		:: bb, b

	CALL RANDOM_SEED()
	CALL RANDOM_NUMBER(bb)
	b = FLOOR((m)*bb)
END SUBROUTINE get_random_vector

PROGRAM Main
	USE sparse_system_solvers
	USE eigv_module

	INTEGER, PARAMETER		:: m = 4	
	REAL, DIMENSION(m,m)	:: A
	REAL, DIMENSION(m)		:: b, x1, x2, x3, lambdar, lambdai
	REAL, DIMENSION(m)		:: rand_vect_real
	REAL, DIMENSION(2*m)	:: rand_idx
	INTEGER i,n
	!REAL start, finish

	! RANDOM MATRIX
	CALL RANDOM_SEED()
	CALL RANDOM_NUMBER(rand_idx)

	CALL RANDOM_SEED()
	CALL RANDOM_NUMBER(rand_vect_real)

	!CALL get_sparse_matrix_llapack(A,m)
	!CALL get_sparse_matrix(A, m, rand_idx, rand_vect_real)
	CALL get_symmetric_sparse_matrix(A, m, rand_idx, rand_vect_real)
	!CALL get_random_vector(b, m)

	! EXAMPLE POSITIVE DEFINITE MATRIX:
	!A = reshape((/2, -2, -3, -2, 5, 4, -3, 4, 5/), (/3,3/))
	! POSITIVE DEFINITE AND SYMMETRIC!
	A = reshape((/2, 1, 0, 0, 1, 3, 2, 0, 0, 2, 5, 2, 0, 0, 2, 3/), (/4,4/))
	!b = (/7,-12,-12/)
	b = (/21, 69, 34, 22/)
	PRINT *, '---------+++---+++-----------'
	CALL solver_v3(A, b, m, x3)
	PRINT *, '---------+++---+++-----------'

	!CALL CPU_TIME(start)	
	!CALL solver_v1(A, b, m, x1)
	!CALL CPU_TIME(finish)
	!WRITE(*, '(" CPU time: ", f12.10, " sec.")') (finish - start)
	
	!CALL CPU_TIME(start)	
	!CALL solver_v2(A, b, m, x2)
	!CALL CPU_TIME(finish)
	!WRITE(*, '(" CPU time: ", f12.10, " sec.")') (finish - start)

	!Compute the eigenvalues to see if matrix is positive definite
	CALL eig(A, lambdar, lambdai)
	PRINT *, '--------------------------'
	PRINT *, 'Real eigenvalues: ', lambdar
	PRINT *, 'Img eigenvalues: ', lambdai
	!PRINT *, '--------------------------'
	!PRINT *, 'LU SOLUTION: ', x1
	!PRINT *, 'CHOLESKY SOLUTION: ', x2

	! COUT matrix to file
	!DO i=1,m
	!	WRITE(12,*) A(i,:)
	!END DO

	!DO i=1,m
	!	WRITE(11,*) b(i)
	!END DO

	!DO i=1,m
	!	WRITE(10,*) x1(i)
	!END DO

	!DO i=1,m
	!	WRITE(9,*) x2(i)
	!END DO
END PROGRAM
