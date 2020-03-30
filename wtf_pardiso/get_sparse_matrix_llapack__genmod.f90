        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 30 14:40:43 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_SPARSE_MATRIX_LLAPACK__genmod
          INTERFACE 
            SUBROUTINE GET_SPARSE_MATRIX_LLAPACK(A,M)
              INTEGER(KIND=4), INTENT(IN) :: M
              REAL(KIND=8) :: A(M,M)
            END SUBROUTINE GET_SPARSE_MATRIX_LLAPACK
          END INTERFACE 
        END MODULE GET_SPARSE_MATRIX_LLAPACK__genmod
