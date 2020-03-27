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
!  Example program for solving symmetric positive definite system of equations.
!  Full case: full functionality of RCI (P)CG is used.
!---------------------------------------------------------------------------
      PROGRAM rci_pcg_f_test_3
        USE MKL_SPBLAS
        IMPLICIT NONE
        INCLUDE 'mkl_rci.fi'
!---------------------------------------------------------------------------
! Define arrays for the upper triangle of the coefficient matrix and
! preconditioner as well as an array for rhs vector
! Compressed sparse row storage is used for sparse representation
!---------------------------------------------------------------------------
        INTEGER N, RCI_request, itercount, expected_itercount, i, info
        PARAMETER (N=8)
        PARAMETER (expected_itercount=8)
        DOUBLE PRECISION  rhs(N)
        INTEGER IA(9)
        INTEGER JA(18)
        DOUBLE PRECISION A(18)
! Fill all arrays containing matrix data.
        DATA IA /1,5,8,10,12,15,17,18,19/
        DATA JA &
&       /1,  3,    6,7,   &
&          2,3,  5,       &
&            3,        8, &
&              4,    7,   &
&                5,6,7,   &
&                  6,  8, &
&                    7,   &
&                      8/
        DATA A &
&       /7.D0,       1.D0,             2.D0, 7.D0,       &
&             -4.D0, 8.D0,       2.D0,                  &
&                    1.D0,                         5.D0,&
&                          7.D0,             9.D0,      &
&                                5.D0, 1.D0, 5.D0,      &
&                                     -1.D0,       5.D0,&
&                                           11.D0,      &
&                                                  5.D0/

        !DOUBLE PRECISION A(4,4)
        !INTEGER IA(5)
        


        !A = reshape((/2, 1, 0, 0, 1, 3, 2, 0, 0, 2, 5, 2, 0, 0, 2, 3/), (/4,4/))
        !b = (/21, 69, 34, 22/)


!---------------------------------------------------------------------------
! Allocate storage for the solver ?par and the initial solution vector
!---------------------------------------------------------------------------
        INTEGER length
        PARAMETER (length=128)
        INTEGER ipar(length)
        DOUBLE PRECISION dpar(length),TMP(N,4)
!---------------------------------------------------------------------------
! Some additional variables to use with the RCI (P)CG solver
!---------------------------------------------------------------------------
        DOUBLE PRECISION  solution(N)
        DOUBLE PRECISION DNRM2, Euclidean_norm, temp(N)
        EXTERNAL DNRM2
        DOUBLE PRECISION alpha, beta
!   Matrix descriptor
        TYPE(MATRIX_DESCR) descrA, descrL
!   CSR matrix representation
        TYPE(SPARSE_MATRIX_T) csrA
        
        PRINT *, A
        PRINT *, "IA", IA
        PRINT *, "JA", JA
        

        !   Create matrix descriptor
        descrA % TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC
        descrA % MODE = SPARSE_FILL_MODE_UPPER
        descrA % DIAG = SPARSE_DIAG_NON_UNIT
        descrL % TYPE = SPARSE_MATRIX_TYPE_DIAGONAL
        descrL % MODE = SPARSE_FILL_MODE_LOWER
        descrL % DIAG = SPARSE_DIAG_NON_UNIT
        alpha = 1.0
        beta  = 0.0
        info = MKL_SPARSE_D_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ONE,N,N,IA,IA(2),JA,A)
        DO I = 1, N
            rhs(I)      = I
            temp(I)     = 0.D0
        END DO

!---------------------------------------------------------------------------
! Initialize the right hand side through matrix-vector product
!---------------------------------------------------------------------------
        !expected_sol = rhs
        !info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,expected_sol,beta,rhs)

        !WRITE(*, *) ' The following r was solved '
        !WRITE(*,800) (rhs(i), i=1,N)
!---------------------------------------------------------------------------
! Initialize the initial guess
!---------------------------------------------------------------------------
        DO I = 1, N
            solution(I)=0.D0
        END DO
!---------------------------------------------------------------------------
! Initialize the solver
!---------------------------------------------------------------------------
        CALL dcg_init(N, solution,rhs, RCI_request,ipar,dpar,TMP)
        IF (RCI_request .NE. 0 ) GO TO 999
!---------------------------------------------------------------------------
! Set the desired parameters:
! INTEGER parameters:
! set the maximal number of iterations to 100
! LOGICAL parameters:
! run the Preconditioned version of RCI (P)CG with preconditioner C_inverse
! DOUBLE PRECISION parameters
! -
!---------------------------------------------------------------------------
        ipar(5)=100
        ipar(11)=1
!---------------------------------------------------------------------------
! Check the correctness and consistency of the newly set parameters
!---------------------------------------------------------------------------
        CALL dcg_check(N,solution,rhs,RCI_request,ipar,dpar,TMP)
        IF (RCI_request .NE. 0 ) GO TO 999
!---------------------------------------------------------------------------
! Compute the solution by RCI (P)CG solver
! Reverse Communications starts here
!---------------------------------------------------------------------------
1       CALL dcg(N,solution,rhs,RCI_request,ipar,dpar,TMP)
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
            CALL DAXPY(N,-1.D0,rhs,1,temp,1)
            Euclidean_norm = DNRM2(N,temp,1)
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
700     CALL dcg_get(N,solution,rhs,RCI_request,ipar,dpar,TMP, itercount)
!---------------------------------------------------------------------------
! Print solution vector: solution(N) and number of iterations: itercount
!---------------------------------------------------------------------------
        WRITE(*, *) ' The system has been solved '
        WRITE(*, *) ' The following solution obtained '
        WRITE(*,800) (solution(i), i=1,N)
800     FORMAT(4(F10.3))
        WRITE(*,900)(itercount)
        WRITE(*, *) ' The following RHS was solved '
        WRITE(*,800) (rhs(i), i=1,N)
        WRITE(*, *) ' The following A was solved '
        WRITE(*,800) (A(i), i=1,N)

900     FORMAT(' Number of iterations: ',1(I2))
        WRITE(*, *) ' The following solution obtained '
        PRINT *, solution
    

!---------------------------------------------------------------------------
! Release internal Intel(R) MKL memory that might be used for computations
! NOTE: It is important to call the routine below to avoid memory leaks
! unless you disable Intel(R) MKL Memory Manager
!---------------------------------------------------------------------------
        CALL MKL_FREE_BUFFERS

        IF (itercount .EQ. expected_itercount .AND. Euclidean_norm .LE. 1.0D-12) THEN
            WRITE( *,'(A,A)') 'This example has successfully PASSED', &
     &  ' through all steps of computation!'
            STOP 0
        ELSE
            WRITE( *,'(A,A,A,I5,A,A,A,E12.5,A)') 'This example may',&
     &  ' have FAILED as either the number of iterations differs from',&
     &  ' the expected number of iterations ',expected_itercount,' or',&
     &  ' the computed solution differs much from the expected',&
     &  ' solution (Euclidean norm is ',Euclidean_norm,'), or both.'
            STOP 1
        END IF
!---------------------------------------------------------------------------
! Release internal Intel(R) MKL memory that might be used for computations
! NOTE: It is important to call the routine below to avoid memory leaks
! unless you disable Intel(R) MKL Memory Manager
!---------------------------------------------------------------------------
999     WRITE( *,'(A,A)') 'This example FAILED as the solver has',&
     &  ' returned the ERROR code', RCI_request
        info = MKL_SPARSE_DESTROY(csrA)
        CALL MKL_FREE_BUFFERS
        STOP 1

      END PROGRAM rci_pcg_f_test_3
