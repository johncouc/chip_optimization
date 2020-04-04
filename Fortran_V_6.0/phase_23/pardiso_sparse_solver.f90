! INSTRUCTIONS TO COMPILE AND LINK ON MAC:
! ifort -i8 -c eigv_module.f90 cg_jacobian_solver.f90 pardiso_sparse_solver.f90 main.f90 -mkl
! ifort -o main eigv_module.o cg_jacobian_solver.o pardiso_sparse_solver.o main.o  -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
!----------------------------------------------------------------------
! Based on "pardiso_sym_f90.f90" example program to show the use of the "PARDISO" routine
! for symmetric linear systems
!---------------------------------------------------------------------
    INCLUDE 'mkl_pardiso.f90'
MODULE pardiso_sparse_solver
    USE mkl_pardiso

    INTEGER(8), PARAMETER :: dp = KIND(1.0D0)
    TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE           :: pt(:)
    integer(8)                            :: iparm(64)
    REAL(KIND=DP):: ddum(1)
    integer(8):: maxfct, mnum,  phase, nrhs, error, msglvl,mtype
    integer(8):: error1, CALLS=0
    integer(8):: i, idum(1)

CONTAINS

   SUBROUTINE pardiso_sym_solver2(acsr,ia,ja, b, m, x)
        integer(8)                                         :: m
        REAL(KIND=DP), DIMENSION(m*m)                     :: b
        REAL(KIND=DP), DIMENSION(m*m),intent(inout)                     :: x

        integer(8), DIMENSION(m*m+1)                         :: ia
        integer(8), DIMENSION(3*m**2-2*m)                         :: ja
        REAL(KIND=DP), DIMENSION(3*m**2-2*m), INTENT(inout)      ::acsr

        !.. Fill all arrays containing matrix data.
        nrhs = 1
        maxfct = 1
        mnum = 1
       ! write(*,*) "PARDISO START"

        !..
        !.. Set up PARDISO control parameter
        !..
        if (CALLS .eq. 0) then
        iparm = 0
        iparm(1) = 1 ! no solver default
        iparm(2) = 2 ! fill-in reordering from METIS
        iparm(4) = 0 !  iterative-direct algorithm - 0 or 62(?)
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
        CALLS=1


        phase=11
        CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m*m, acsr, ia, ja,idum, nrhs, iparm, msglvl, ddum, ddum, error)

                end if

     
        !.. Reordering and Symbolic Factorization, This step also allocates
        ! all memory that is necessary for the factorization
        phase = 11 ! only reordering and symbolic factorization
        
        IF (error /= 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', error
        END IF


        !change between iterative-direct here:


        if (CALLS .gt. 2 ) then
        IPARM(4)=0
        end if
        
        phase = 23 !  factorization & solve
        CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m*m, acsr, ia, ja, &
                      idum, nrhs, iparm, msglvl, b, x, error)
        IF (error /= 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', error
        ENDIF


        !check cg iterative here:
       ! if (iparm(20) .lt. 1) then
       !         write (*,*)     "PARDISO CG ITERS:", iparm(20)
       ! endif



        IF (error /= 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', error
       !    GOTO 1000
        ENDIF
        1000 CONTINUE
        END SUBROUTINE pardiso_sym_solver2

        SUBROUTINE pardiso_sym_solver3(acsr,ia,ja, b, m, x)
        integer(8)                                         :: m 
        REAL(KIND=DP), DIMENSION(m*m)                     :: b
        REAL(KIND=DP), DIMENSION(m*m),intent(inout)                     :: x

        !   TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE           :: pt(:)
        !.. All other variables
        integer(8), DIMENSION(m*m+1)                         :: ia
        integer(8), DIMENSION(3*m**2-2*m)                         :: ja
        REAL(KIND=DP), DIMENSION(3*m**2-2*m), INTENT(inout)      ::acsr

       !        integer(8)::mtype
        !.. Fill all arrays containing matrix data.
        nrhs = 1
        maxfct = 1
        mnum = 1

        phase = 33 ! only solving
        CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m*m, acsr, ia, ja, &
                      idum, nrhs, iparm, msglvl, b, x, error)

        CALLS=CALLS+1
!        write(*,*) "end of pardiso, iparm(20) =", iparm(20)
!        read(*,*)        
        END SUBROUTINE pardiso_sym_solver3
        
        SUBROUTINE pardiso_release(m)
        integer(8) :: m
           
          write (*,*) "IPARM(20) : ", iparm(20), " IPARM(4) : ",iparm(4)
           
                phase = -1 ! release internal memory
        CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m*m, ddum, idum, idum,idum, nrhs, iparm, msglvl, ddum, ddum, error)
         
        !.. Termination and release of memory
       
       ! IF (ALLOCATED(jaa)) DEALLOCATE(jaa)
       ! IF (ALLOCATED(acsrr)) DEALLOCATE(acsrr)
       ! IF (ALLOCATED(iparm)) DEALLOCATE(iparm)
        if (ALLOCATED(pt)) DEALLOCATE(pt)

        IF (error1 /= 0) THEN
           WRITE(*,*) 'The following ERROR on release stage was detected: ', error
           STOP 1
        ENDIF

        IF (error /= 0) STOP 1
    END SUBROUTINE pardiso_release















END MODULE pardiso_sparse_solver
