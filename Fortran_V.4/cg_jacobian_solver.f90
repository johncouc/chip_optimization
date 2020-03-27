!===============================================================================
! Copyright 2005-2019 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

!  Content: Intel(R) MKL RCI (P)CG Fortran example
!
!*******************************************************************************

!---------------------------------------------------------------------------
!  Based on "cg_jacobi_precon.f90" example program for solving
!               symmetric positive definite system of equations.
!  Full case: full functionality of RCI (P)CG is used.
!---------------------------------------------------------------------------

MODULE cg_jacobian_solver
    USE MKL_SPBLAS
    INCLUDE 'mkl_rci.fi'
CONTAINS
    SUBROUTINE  uncompressed_to_CSR_converter(ADNS, m, JA, IA, Acsr)
    ! Convert a sparse matrix in uncompressed representation to the CSR format
    ! Double precision
        IMPLICIT NONE
        INTEGER                 :: m, info
        REAL(8), DIMENSION(m,m) :: ADNS
        REAL(8), DIMENSION(m*m) :: Acsr
        INTEGER, DIMENSION(m+1) :: IA
        INTEGER, DIMENSION(m*m) :: JA
        INTEGER, DIMENSION(8)   :: job

        info = 0
        job(1)=0
        job(2)= 1
        job(3)= 1
        job(4)= 1
        job(5)=m*m
        job(6)=1

        Acsr = 0
        JA = 0
        IA = 0

        CALL mkl_ddnscsr(job,m,m,ADNS,m,Acsr,JA,IA,info)
    END SUBROUTINE uncompressed_to_CSR_converter

    SUBROUTINE pardiso_solver(A, b, m, solution)

            IMPLICIT NONE
            INTEGER                                         :: m, zero_index
            REAL(8), DIMENSION(m,m), INTENT(inout)          :: A
            REAL(8), DIMENSION(m), INTENT(in)               :: b
            REAL(8), DIMENSION(m)                           :: solution
            REAL(8), DIMENSION(m*m)                         :: acsr
            INTEGER, DIMENSION(m+1)                         :: IA
            INTEGER, DIMENSION(m*m)                         :: JA
            INTEGER, DIMENSION(:), ALLOCATABLE              :: JAA
            REAL(8), DIMENSION(:), ALLOCATABLE              :: acsrr

            INTEGER RCI_request, itercount, i, info
    !---------------------------------------------------------------------------
    ! Allocate storage for the solver ?par and the initial solution vector
    !---------------------------------------------------------------------------
            INTEGER length
            PARAMETER (length=128)
            INTEGER ipar(length)
            REAL(8) dpar(length),TMP(m,4)
    !---------------------------------------------------------------------------
    ! Some additional variables to use with the RCI (P)CG solver
    !---------------------------------------------------------------------------
            REAL(8) DNRM2, Euclidean_norm, temp(m)
            EXTERNAL DNRM2
            REAL(8) alpha, beta
            ! Matrix descriptor
            TYPE(MATRIX_DESCR) descrA, descrL
            ! CSR matrix representation
            TYPE(SPARSE_MATRIX_T) csrA
            ! Create matrix descriptor
            descrA % TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC
            descrA % MODE = SPARSE_FILL_MODE_UPPER
            descrA % DIAG = SPARSE_DIAG_NON_UNIT
            descrL % TYPE = SPARSE_MATRIX_TYPE_DIAGONAL
            descrL % MODE = SPARSE_FILL_MODE_LOWER
            descrL % DIAG = SPARSE_DIAG_NON_UNIT

            ! Convert a sparse matrix in uncompressed representation to the CSR format
            CALL uncompressed_to_CSR_converter(A, m, JA, IA, acsr)

            ! Remove zeros
            ! Get the last element of the array and substract 1 as fortran is 0-based.
            zero_index = IA(size(IA)) - 1
            ALLOCATE(JAA(zero_index))
            ALLOCATE(acsrr(zero_index))

            JAA = JA(:zero_index)
            acsrr = acsr(:zero_index)

            alpha = 1.0
            beta  = 0.0
            info = MKL_SPARSE_D_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ONE,m,m,IA,IA(2),JAA,acsrr)
            temp     = 0.D0

    !---------------------------------------------------------------------------
    ! Initialize the right hand side through matrix-vector product
    !---------------------------------------------------------------------------
    !       info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,expected_sol,beta,rhs)
    !---------------------------------------------------------------------------
    ! Initialize the initial guess
    !---------------------------------------------------------------------------
            DO I = 1, m
                solution(I)=0.D0
            END DO
    !---------------------------------------------------------------------------
    ! Initialize the solver
    !---------------------------------------------------------------------------
            CALL dcg_init(m, solution,b, RCI_request,ipar,dpar,TMP)
            IF (RCI_request .NE. 0 ) GO TO 999
    !---------------------------------------------------------------------------
    ! Set the desired parameters:
    ! INTEGER parameters:
    ! set the maximal number of iterations to 100
    ! LOGICAL parameters:
    ! run the Preconditioned version of RCI (P)CG with preconditioner C_inverse
    ! DOUBLE PRECISION parameters
    !---------------------------------------------------------------------------
            ipar(5)=100
            ipar(11)=1
    !---------------------------------------------------------------------------
    ! Check the correctness and consistency of the newly set parameters
    !---------------------------------------------------------------------------
            CALL dcg_check(m,solution,b,RCI_request,ipar,dpar,TMP)
            IF (RCI_request .NE. 0 ) GO TO 999
    !---------------------------------------------------------------------------
    ! Compute the solution by RCI (P)CG solver
    ! Reverse Communications starts here
    !---------------------------------------------------------------------------
1         CALL dcg(m,solution,b,RCI_request,ipar,dpar,TMP)
    !---------------------------------------------------------------------------
    ! If RCI_request=0, then the solution was found according to the requested
    ! stopping tests. In this case, this means that it was found after 100
    ! iterations.
    !---------------------------------------------------------------------------
                IF (RCI_request .EQ. 0) THEN
                GO TO 700
    !---------------------------------------------------------------------------
    ! If RCI_request=1, then compute the vector A*TMP(:,1)
    ! and put the result in vector TMP(:,2)
    !---------------------------------------------------------------------------
            ELSE IF (RCI_request .EQ. 1) THEN
            info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,TMP,beta,TMP(1,2))
                GO TO 1
    !---------------------------------------------------------------------------
    ! If RCI_request=2, then do the user-defined stopping test: compute the
    ! Euclidean norm of the actual residual using Intel(R) MKL routines and check if
    ! it is less than 1.D-8
    !---------------------------------------------------------------------------
            ELSE IF (RCI_request .EQ. 2) THEN
                info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,solution,beta,temp)
                CALL DAXPY(m,-1.D0,b,1,temp,1)
                Euclidean_norm = DNRM2(m,temp,1)
                IF (Euclidean_norm .GT. 1.D-8) THEN
    !---------------------------------------------------------------------------
    ! The solution has not been found yet according to the user-defined stopping
    ! test. Continue RCI (P)CG iterations.
    !---------------------------------------------------------------------------
                    GO TO 1
                ELSE
    !---------------------------------------------------------------------------
    ! The solution has been found according to the user-defined stopping test
    !---------------------------------------------------------------------------
                    GO TO 700
                END IF
    !---------------------------------------------------------------------------
    ! If RCI_request=3, then compute apply the preconditioner matrix C_inverse
    ! on vector TMP(:,3) and put the result in vector TMP(:,4)
    !---------------------------------------------------------------------------
            ELSE IF (RCI_request .EQ. 3) THEN
                 info = MKL_SPARSE_D_TRSV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,&
    &                                     descrL,TMP(1,3),TMP(1,4))
                GO TO 1
            ELSE
    !---------------------------------------------------------------------------
    ! If RCI_request=anything else, then dcg subroutine failed
    ! to compute the solution vector: solution(N)
    !---------------------------------------------------------------------------
                GO TO 999
            END IF
    !---------------------------------------------------------------------------
    ! Reverse Communication ends here
    ! Get the current iteration number
    !---------------------------------------------------------------------------
700         CALL dcg_get(m,solution,b,RCI_request,ipar,dpar,TMP, itercount)
    !---------------------------------------------------------------------------
    ! Release internal Intel(R) MKL memory that might be used for computations
    ! NOTE: It is important to call the routine below to avoid memory leaks
    ! unless you disable Intel(R) MKL Memory Manager
    !---------------------------------------------------------------------------
            CALL MKL_FREE_BUFFERS
999         info = MKL_SPARSE_DESTROY(csrA)
            CALL MKL_FREE_BUFFERS

            DEALLOCATE(JAA)
            DEALLOCATE(acsrr)
    END SUBROUTINE pardiso_solver
END MODULE cg_jacobian_solver
