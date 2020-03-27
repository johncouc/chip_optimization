! COMMANDS TO RUN:
! gfortran -c main.f90 eigv_module.f90 fvm.f90
! gfortran -o main.f90 eigv_module.f90 fvm.f90
! ./main
! The matrices and results will be saved in output files 'fort.*'
! where * respresents 11,12 for the input matrix and vector
! 9 and 10 - output solutions for two methods

! Module is based on: https://www.netlib.org/lapack/lug/node38.html
MODULE sparse_system_solvers
	IMPLICIT NONE

	SAVE
        public solver_v1,solver_v11
        integer, dimension(14400) ::IPIVV
CONTAINS

	! Solves general band matrices including tridiagonal matrices:
	! LU factorization with partial pivoting
	SUBROUTINE solver_v1(A, b, m, x)
		INTEGER :: m
		REAL(8), DIMENSION(m,m), INTENT(inout)		:: A
	!	REAL(8), DIMENSION(m,m)					:: AA,BB
		REAL(8), DIMENSION(m), INTENT(in)			:: b
		REAL(8), DIMENSION(m), INTENT(out)			:: x
		INTEGER, DIMENSION(m)					:: IPIV
		INTEGER info, l_work,k
	!	write (*,*) "SPARSE:::: ) " , maxval(A)
	!	write (*,*) "SPARSE:::: ) " , b
		! Make matrix and vector copies
        !	AA = A
		x = b

		!Compute LU factorization 
        !CALL DGETRF(SIZE(AA,1), SIZE(AA,2), AA, m, IPIV, info)
        CALL DGETRF(SIZE(A,1), SIZE(A,2), A, m, IPIV, info)
        do k=1,m
        IPIVV(k)=IPIV(k)
        enddo
         

		! Solve system
		CALL DGETRS('N',SIZE(A, DIM=2),1, A, SIZE(A, DIM=1), IPIV, x, SIZE(x,1), info)
        !        write(*,*) x
        !        read(*,*)         


	END SUBROUTINE solver_v1
	SUBROUTINE solver_v11(A, b, m, x)
		INTEGER								  	:: m
		REAL(8), DIMENSION(m,m), INTENT(in)		:: A
		REAL(8), DIMENSION(m), INTENT(in)			:: b
		REAL(8), DIMENSION(m), INTENT(out)			:: x
		INTEGER, DIMENSION(m)					:: IPIV
		INTEGER info, l_work,k
		
                ! Make matrix and vector copies
        	x = b

                do k=1,m
                IPIV(k)=IPIVV(k)
                enddo
                
               ! read(*,*)
		CALL DGETRS('N',SIZE(A, DIM=2),1, A, SIZE(A, DIM=1), IPIV, x, SIZE(x,1), info)
               ! write(*,*) x
               ! read(*,*) 
	END SUBROUTINE solver_v11

	! Solves symmetric and positive definite matrices including band matrices:
	! Cholesky factorization
	SUBROUTINE solver_v2(A, b, m, x)
		INTEGER								  	:: m
		REAL(8), DIMENSION(m,m), INTENT(inout)		:: A
		!REAL(8), DIMENSION(m,m)					:: AA, A_u
		REAL(8), DIMENSION(m), INTENT(in)			:: b
		REAL(8), DIMENSION(m), INTENT(inout)			:: x
!		INTEGER, DIMENSION(m)					:: IPIV
		INTEGER info, l_work
		
		! Make matrix and vector copies
!		AA = A
		x = b

		! Solves a system of linear equations A*X = B where 
		! A is an N-by-N symmetric positive definite matrix
		! Performs Cholesky factorization
		CALL DPOTRF('L', SIZE(A,1), A, SIZE(A,1), info)
	        CALL DPOTRS('L', SIZE(A,1), 1, A, SIZE(A,1), x, SIZE(x,1), info)
!                write(*,*) x
!                read(*,*)
	END SUBROUTINE solver_v2

	SUBROUTINE solver_v22(A, b, m, x)
		INTEGER								  	:: m
		REAL(8), DIMENSION(m,m), INTENT(in)		:: A
		!REAL(8), DIMENSION(m,m)					:: AA, A_u
		REAL(8), DIMENSION(m), INTENT(in)			:: b
		REAL(8), DIMENSION(m), INTENT(inout)			:: x
!		INTEGER, DIMENSION(m)					:: IPIV
		INTEGER info, l_work
		
		! Make matrix and vector copies
!		AA = A
		x = b

		! Solves a system of linear equations A*X = B where 
		! A is an N-by-N symmetric positive definite matrix
		! Performs Cholesky factorization
	!	CALL DPOTRF('L', SIZE(A,1), 1, A, SIZE(A,1), x, SIZE(x,1), info)
	        CALL DPOTRS('L', SIZE(A,1), 1, A, SIZE(A,1), x, SIZE(x,1), info)
!                write(*,*) x
!                read(*,*)
	END SUBROUTINE solver_v22


	! BASED ON THIS COMMENT: https://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=5048
	! TRY BAND SOLVERS IN LLAPACK
	! Computes the solution to a real system of linear equations A * X = B,
 	! where A is an N-by-N symmetric positive definite band matrix and X
 	! and B are N-by-NRHS matrices.
	SUBROUTINE solver_v3(A, b, m, x)
		INTEGER								  	:: m, nsuper, info, i, j
		REAL(8), DIMENSION(m,m), INTENT(in)		:: A
		REAL(8), DIMENSION(m,m)					:: A_u
		REAL(8), DIMENSION(CEILING(m/2.0),m)			:: A_banded
		REAL(8), DIMENSION(m)						:: b, x

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

SUBROUTINE get_sparse_matrix(A, m, rand_idx, rand_vect_real)
	INTEGER, INTENT(in) 				:: m
	REAL(8), DIMENSION(m,m)				:: A
	REAL(8), DIMENSION(m), INTENT(IN)		:: rand_vect_real
	REAL(8), DIMENSION(2*m), INTENT(IN)	:: rand_idx
	REAL(8), DIMENSION(2*m)				:: rand_idx_COPY
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
	REAL(8), DIMENSION(m,m)				:: A
	REAL(8), DIMENSION(m)					:: rand_vect_real
	REAL(8), DIMENSION(2*m), INTENT(IN)	:: rand_idx
	REAL(8), DIMENSION(2*m)				:: rand_idx_COPY
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
	REAL(8), DIMENSION(m)		:: bb, b

	CALL RANDOM_SEED()
	CALL RANDOM_NUMBER(bb)
	b = FLOOR((m)*bb)
END SUBROUTINE get_random_vector



