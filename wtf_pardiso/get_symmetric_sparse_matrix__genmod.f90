        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 30 14:40:43 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_SYMMETRIC_SPARSE_MATRIX__genmod
          INTERFACE 
            SUBROUTINE GET_SYMMETRIC_SPARSE_MATRIX(A,M,RAND_IDX,        &
     &RAND_VECT_REAL)
              INTEGER(KIND=4), INTENT(IN) :: M
              REAL(KIND=8) :: A(M,M)
              REAL(KIND=8), INTENT(IN) :: RAND_IDX(2*M)
              REAL(KIND=8) :: RAND_VECT_REAL(M)
            END SUBROUTINE GET_SYMMETRIC_SPARSE_MATRIX
          END INTERFACE 
        END MODULE GET_SYMMETRIC_SPARSE_MATRIX__genmod
