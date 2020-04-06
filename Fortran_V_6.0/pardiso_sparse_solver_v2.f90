INCLUDE 'mkl_pardiso.f90'
MODULE pardiso_sparse_solver_v2
    USE mkl_pardiso

    TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE       :: pt(:)
    INTEGER(8) :: maxfct, mnum, mtype, phase, msglvl, error, nrhs
    INTEGER(8) :: iteration_count = 0, MODE=1 ! set MODE=0 to use direct/ Mode=1 to use iterative
    INTEGER(8), dimension(:), ALLOCATABLE       :: perm
    INTEGER(8), dimension(64)                   :: iparm
    REAL(8), DIMENSION(1)                       :: ddum

CONTAINS
SUBROUTINE pardiso_sym_solver2(acsr, ia, ja, b, m, x)
    INTEGER(8)                              :: m
    REAL(8), DIMENSION(3*m**2-2*m)          :: acsr
    REAL(8), dimension(m*m)                 :: b, x
    INTEGER(8), dimension(3*m**2-2*m)       :: ja
    INTEGER(8), dimension(m*m+1)            :: ia
    INTEGER                                 :: i

     iparm = 0
     iparm(1) = 1
     iparm(2) = 2
     iparm(3) = 1
     iparm(4) = 0
     iparm(5) = 0
     iparm(6) = 0
     iparm(7) = 0
     iparm(8) = 2
     iparm(9) = 0
     iparm(10) = 13
     iparm(11) = 1
     iparm(12) = 0
     iparm(13) = 1
     iparm(14) = -1
     iparm(15) = 0
     iparm(16) = 0
     iparm(17) = -1
     iparm(18) = -1
     iparm(19) = -1
     iparm(20) = -1
     iparm(28) = 0! single precision

    IF(iteration_count == 0) THEN
        ALLOCATE(perm(m)); perm=0
        ALLOCATE (pt(64))
        DO i=1,64
          pt(i)%DUMMY =  0
        ENDDO

        maxfct= 1
        mnum = 1
        mtype = 2
        msglvl= 0
        nrhs = 1
    END IF

    IF(iteration_count == 0 .OR. MODE == 0) THEN
        phase=11
        CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m*m, acsr, ia, &
                     ja, perm, nrhs, iparm, msglvl, b, x, error)

        phase=23
        CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m*m, acsr, ia, &
                     ja, perm, nrhs, iparm, msglvl, b, x, error)
        iteration_count = iteration_count + 1

    ELSE
        phase=23
        iparm(4)=62
        CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m*m, acsr, ia, &
                     ja, perm, nrhs, iparm, msglvl, b, x, error)
    ENDIF

    !write(*,'(a)') "Solve completed ..."
    !write (*,*) "IPARM(20) : ", iparm(20), " IPARM(4) : ",iparm(4)

END SUBROUTINE

SUBROUTINE pardiso_sym_solver3(acsr, ia, ja, b, m, x)
    INTEGER(8)                              :: m
    REAL(8), DIMENSION(3*m**2-2*m)          :: acsr
    REAL(8), dimension(m*m)                 :: b, x
    INTEGER(8), dimension(3*m**2-2*m)       :: ja
    INTEGER(8), dimension(m*m+1)            :: ia

    phase = 33 ! only solving
    CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m*m, acsr, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)

END SUBROUTINE


SUBROUTINE pardiso_release(m)
    INTEGER(8) :: m
                      
    phase = -1 ! release internal memory
    CALL pardiso_64(pt, maxfct, mnum, mtype, phase, m*m, ddum, perm, perm, perm, nrhs, iparm, msglvl, ddum, ddum, error)
     
    !.. Termination and release of memory
    IF (ALLOCATED(pt)) DEALLOCATE(pt)
    IF (ALLOCATED(perm)) DEALLOCATE(perm)


    IF (error1 /= 0) THEN
       WRITE(*,*) 'The following ERROR on release stage was detected: ', error
       STOP 1
    ENDIF

    IF (error /= 0) STOP 1
END SUBROUTINE pardiso_release
END MODULE
