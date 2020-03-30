! INSTRUCTIONS TO COMPILE AND LINK ON MAC:
! ifort -i8 -c eigv_module.f90 cg_jacobian_solver.f90 pardiso_sparse_solver.f90 main.f90 -mkl
! ifort -o main eigv_module.o cg_jacobian_solver.o pardiso_sparse_solver.o main.o  -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
!----------------------------------------------------------------------
! Based on "pardiso_sym_f90.f90" example program to show the use of the "PARDISO" routine
! for symmetric linear systems
!---------------------------------------------------------------------
MODULE pardiso_sparse_solver
    USE mkl_pardiso
    !INCLUDE 'mkl_pardiso.fi'
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
CONTAINS
    SUBROUTINE  uncompressed_to_CSR_converter(ADNS, m, ja, ia, Acsr)
    ! Convert a sparse matrix in uncompressed representation to the CSR format
    ! Double precision
        IMPLICIT NONE
        INTEGER                 :: m, info
        REAL(8), DIMENSION(m,m) :: ADNS
        REAL(8), DIMENSION(m*m) :: Acsr
        INTEGER, DIMENSION(m+1) :: ia
        INTEGER, DIMENSION(m*m) :: ja
        INTEGER, DIMENSION(8)   :: job

        info = 0
        job(1)=0
        job(2)= 1
        job(3)= 1
        job(4)= 1
        job(5)=m*m
        job(6)=1

        Acsr = 0
        ja = 0
        ia = 0

        CALL mkl_ddnscsr(job,m,m,ADNS,m,Acsr,ja,ia,info)
    END SUBROUTINE uncompressed_to_CSR_converter
    
    SUBROUTINE pardiso_sym_solver(A, b, m, x)
        INTEGER                                         :: m, zero_index
        REAL(KIND=DP), DIMENSION(m,m), INTENT(inout)    :: A
        REAL(KIND=DP), DIMENSION(m)                     :: b
        REAL(KIND=DP), DIMENSION(m)                     :: x
        !.. Internal solver memory pointer
        TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE           :: pt(:)
        !.. All other variables
        INTEGER, ALLOCATABLE                            :: iparm(:)
        REAL(KIND=DP), DIMENSION(m*m)                   :: acsr
        INTEGER, DIMENSION(m+1)                         :: ia
        INTEGER, DIMENSION(m*m)                         :: ja
        INTEGER, DIMENSION(:), ALLOCATABLE              :: JAA
        REAL(KIND=DP), DIMENSION(:), ALLOCATABLE        :: acsrr
        REAL(KIND=DP) ddum(1)
        INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
        INTEGER error1
        INTEGER i, idum(1)

        !.. Fill all arrays containing matrix data.
        nrhs = 1
        maxfct = 1
        mnum = 1

        ! Convert a sparse matrix in uncompressed representation to the CSR format
        CALL uncompressed_to_CSR_converter(A, m, ja, ia, acsr)

        zero_index = ia(size(ia)) - 1
        ALLOCATE(jaa(zero_index))
        ALLOCATE(acsrr(zero_index))
        jaa = ja(:zero_index)
        acsrr = acsr(:zero_index)

        !..
        !.. Set up PARDISO control parameter
        !..
        ALLOCATE(iparm(64))

        iparm = 0
        iparm(1) = 1 ! no solver default
        iparm(2) = 2 ! fill-in reordering from METIS
        iparm(4) = 0 ! no iterative-direct algorithm
        iparm(5) = 0 ! no user fill-in reducing permutation
        iparm(6) = 0 ! =0 solution on the first n components of x
        iparm(8) = 2 ! numbers of iterative refinement steps
        iparm(10) = 13 ! perturb the pivot elements with 1E-13
        iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
        iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
        iparm(14) = -1 ! Output: number of perturbed pivots
        iparm(17) = -1 ! Output: number of non-zero elements in the factors
        iparm(18) = -1 ! Output: number of nonzeros in the factor LU
        iparm(19) = -1 ! Output: Mflops for LU factorization
        iparm(20) = -1 ! Output: Numbers of CG Iterations
        iparm(36) = 0 ! matrix stored in CSR format
        

        error  = 0 ! initialize error flag
        msglvl = 0 ! print statistical information
        mtype  = 2 ! symmetric, positive definite

        !.. Initialize the internal solver memory pointer. This is only
        ! necessary for the FIRST call of the PARDISO solver.

        ALLOCATE (pt(64))
        DO i = 1, 64
           pt(i)%DUMMY =  0
        END DO

        !.. Reordering and Symbolic Factorization, This step also allocates
        ! all memory that is necessary for the factorization

        phase = 11 ! only reordering and symbolic factorization
        CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m, acsrr, ia, jaa, &
                      idum, nrhs, iparm, msglvl, ddum, ddum, error)
        
        !WRITE(*,*) 'Reordering completed ... '
        IF (error /= 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', error
           GOTO 1000
        END IF
        !WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
        !WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)

        !.. Factorization.
        phase = 22 ! only factorization
        CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m, acsrr, ia, jaa, &
                      idum, nrhs, iparm, msglvl, ddum, ddum, error)
        !WRITE(*,*) 'Factorization completed ... '
        IF (error /= 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', error
           GOTO 1000
        ENDIF

        !.. Back substitution and iterative refinement
        iparm(8) = 2 ! max numbers of iterative refinement steps
        phase = 33 ! only solving
        CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m, acsrr, ia, jaa, &
                      idum, nrhs, iparm, msglvl, b, x, error)
        !WRITE(*,*) 'Solve completed ... '
        IF (error /= 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', error
           GOTO 1000
        ENDIF
          
        1000 CONTINUE
        !.. Termination and release of memory
        phase = -1 ! release internal memory
        CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m, ddum, idum, idum, &
                      idum, nrhs, iparm, msglvl, ddum, ddum, error1)

        IF (ALLOCATED(jaa)) DEALLOCATE(jaa)
        IF (ALLOCATED(acsrr)) DEALLOCATE(acsrr)
        IF (ALLOCATED(iparm)) DEALLOCATE(iparm)

        IF (error1 /= 0) THEN
           WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
           STOP 1
        ENDIF

        IF (error /= 0) STOP 1
    END SUBROUTINE pardiso_sym_solver
END MODULE pardiso_sparse_solver
