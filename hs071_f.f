C Copyright (C) 2002, 2010 Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Eclipse Public License.
C
C    $Id$
C
C =============================================================================
C
C =============================================================================
C
C                            Main driver program
C
C =============================================================================
C
      PROGRAM example
      USE timings
      USE fvm

      IMPLICIT NONE

C     include the Ipopt return codes
      INCLUDE 'IpReturnCodes.inc'
      
C     Size of the problem (number of variables and equality constraints)
      INTEGER(8) N, M, NELE_JAC, NELE_HESS, IDX_STY
      PARAMETER  (N =100*100, M = 1, NELE_JAC = N, NELE_HESS = 0)
      PARAMETER  (IDX_STY = 1 )

C     Space for multipliers and constraints
      DOUBLE PRECISION LAM(M)
      DOUBLE PRECISION G(M)

C     Vector of variables
      DOUBLE PRECISION X(N)

C     Vector of lower and upper bounds
      DOUBLE PRECISION X_L(N), X_U(N), Z_L(N), Z_U(N)
      DOUBLE PRECISION G_L(M), G_U(M)

C     Private data for evaluation routines
C     This could be used to pass DOUBLE PRECISION and INTEGER(8) arrays untouched
C     to the evaluation subroutines EVAL_*
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE ::DAT
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE:: IDAT

C     Place for storing the Ipopt Problem Handle
C     for 32 bit platforms
C      INTEGER IPROBLEM
C      INTEGER IPCREATE
C     for 64 bit platforms:
      INTEGER*8 IPROBLEM
      INTEGER*8 IPCREATE

      INTEGER(8) IERR
      INTEGER(8) IPSOLVE, IPADDSTROPTION
      INTEGER(8) IPADDNUMOPTION, IPADDINTOPTION
      INTEGER(8) IPOPENOUTPUTFILE

      DOUBLE PRECISION F
      INTEGER(8) i

      REAL(8), DIMENSION(1,N) :: v , GRAD

      CALL mkl_set_num_threads(4)
      CALL startClock()
      CALL tic()

      d=((floor(sqrt(N*1d0))+1))
      
      ALLOCATE (DAT(3* d**2-2*d,2))
      ALLOCATE (IDAT (d+1,d+1))
      ALLOCATE (DAT3(d**2+1))
      ALLOCATE (DAT4(3*d**2-2*d))
      DAT=0
      IDAT=0 
      DAT3=0
      DAT4=0
      v=0.4d0 
 
C     Set initial point and bounds:
      X =  0.4
      X_L =0
      X_U= 1

C     Set bounds for the constraints
       G_L = 0d0
       G_U = (0.399d0)*(0.005d0)**2
      
C     First create a handle for the Ipopt problem (and read the options file)
      IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC,0, 
     1     IDX_STY, EV_F_CHIP, EV_G, EV_GRAD_F_CHIP, EV_JAC_G)
      IF (IPROBLEM.eq.0) THEN
         WRITE(*,*) 'Error creating an Ipopt Problem handle.'
         stop
      endif

C     Open an output file
      IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 5)
      IF (IERR .NE.0 ) THEN
         WRITE(*,*) 'Error opening the Ipopt output file.'
         goto 9000
      endif

C     Set a string option
      IERR = IPADDSTROPTION(IPROBLEM, 'mu_strategy', 'adaptive')
      IF (IERR .NE.0 ) goto 9990
      IERR = IPADDSTROPTION(IPROBLEM,'hessian_approximation',
     1'limited-memory')


      IF (IERR .NE.0 ) goto 9990

C     Set an INTEGER option
      IERR = IPADDINTOPTION(IPROBLEM, 'max_iter',10)
      IF (IERR .NE.0 ) goto 9990

      !IERR = IPADDINTOPTION(IPROBLEM, 'print_level',0)
      !IF (IERR .NE.0 ) goto 9990

C     Set a DOUBLE PRECISION option
      IERR = IPADDNUMOPTION(IPROBLEM, 'tol', 1.d-7)
      IF (IERR .NE.0 ) goto 9990

      IERR = IPADDNUMOPTION(IPROBLEM, 'derivative_test_perturbation',
     1 5.d-5)
       IERR = IPADDNUMOPTION(IPROBLEM, 'derivative_test_tol',
     1 1.d-3)
       IERR = IPADDNUMOPTION(IPROBLEM, 'point_perturbation_radius',
     1 1.d-4) 
        IERR = IPADDNUMOPTION(IPROBLEM,'expect_infeasible_problem_ctol',
     1 1.d-7)
      p=3d0

C     CALL optimization routine
      IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, DAT, IDAT)
      CALL stopClock()
      CALL toc()

 
      READ(*,*)
      OPEN(unit = 2, file = "optimization_result")
      DO I=1,N
        WRITE (2,*) X(I)
      ENDDO
      CLOSE(2)
      DEALLOCATE(DAT,IDAT,DAT3,DAT4)

C     Output:
      IF( IERR .EQ. IP_SOLVE_SUCCEEDED ) THEN
         WRITE(*,*)
         WRITE(*,*) 'The solution was found.'
         WRITE(*,*)
         WRITE(*,*) 'The final value of the objective function is ',F
         WRITE(*,*)
         WRITE(*,*) 'The optimal values of X are:'
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,*) 'The error code is ',IERR
         WRITE(*,*)
      ENDIF

 9000 CONTINUE

C     Clean up
      CALL IPFREE(IPROBLEM)
      CALL pardiso_release(d)
      STOP
C
 9990 CONTINUE
      WRITE(*,*) 'Error setting an option'
      GOTO 9000

      END
