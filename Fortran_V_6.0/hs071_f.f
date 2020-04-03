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
      program example

C
      use timings
      use fvm

      implicit none
C
C     include the Ipopt return codes
C
      include 'IpReturnCodes.inc'
      
C     Size of the problem (number of variables and equality constraints)
C
      integer(8)     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
      parameter  ( N =160*160, M = 1, NELE_JAC = N, NELE_HESS = 0)
      parameter  (IDX_STY = 1 )
      !REAL(8) start, finish

C
C     Space for multipliers and constraints
C
      double precision LAM(M)
      double precision G(M)
C
C     Vector of variables
C
      double precision X(N)
C
C     Vector of lower and upper bounds
C
      double precision X_L(N), X_U(N), Z_L(N), Z_U(N)
      double precision G_L(M), G_U(M)
C
C     Private data for evaluation routines
C     This could be used to pass double precision and integer(8) arrays untouched
C     to the evaluation subroutines EVAL_*
C
      double precision, dimension (:,:), allocatable ::DAT
      double precision, dimension (:,:), allocatable:: IDAT
C     Place for storing the Ipopt Problem Handle
C
CC     for 32 bit platforms
C      integer IPROBLEM
C      integer IPCREATE
C     for 64 bit platforms:
      integer*8 IPROBLEM
      integer*8 IPCREATE
C
      integer(8) IERR
      integer(8) IPSOLVE, IPADDSTROPTION
      integer(8) IPADDNUMOPTION, IPADDINTOPTION
      integer(8) IPOPENOUTPUTFILE
C
      double precision F
C      integer(8) i
C

      real(8), dimension(1,N) :: v , GRAD
            !  call fvm_simulate(v, d,DAT,IDAT)
C     d**2,d**2+1

      call mkl_set_num_threads(4)
      call startClock()
      call tic()

      !CALL tic()
      !CALL CPU_TIME(start)

      d=((floor(sqrt(N*1d0))+1))
      
      allocate (DAT(3* d**2-2*d,2))
      allocate (IDAT (d+1,d+1))
      allocate (DAT3(d**2+1))
      allocate (DAT4(3*d**2-2*d))
      DAT=0
      IDAT=0 
      DAT3=0
      DAT4=0
      v=0.4d0 
 
  

C     Set initial point and bounds:
C
      X =  0.4
C        open (unit = 2, file = "optimization_result")
C      do I=1,N
C      read (2,*) X(I)
C      enddo
C        close(2)
      X_L =0
      X_U= 1
C
C     Set bounds for the constraints
C
       G_L = 0d0
  
       G_U = (0.399d0)*(0.005d0)**2
      
C     First create a handle for the Ipopt problem (and read the options
C     file)
C
!["ipopt.print_level"] = 0

      IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC,0, 
     1     IDX_STY, EV_F_CHIP, EV_G, EV_GRAD_F_CHIP, EV_JAC_G)
      if (IPROBLEM.eq.0) then
         write(*,*) 'Error creating an Ipopt Problem handle.'
         stop
      endif
C
C     Open an output file
C
      IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 5)
      if (IERR.ne.0 ) then
         write(*,*) 'Error opening the Ipopt output file.'
         goto 9000
      endif
C
C     Note: The following options are only examples, they might not be
C           suitable for your optimization problem.
C
C     Set a string option
C
      IERR = IPADDSTROPTION(IPROBLEM, 'mu_strategy', 'adaptive')
      if (IERR.ne.0 ) goto 9990
      IERR = IPADDSTROPTION(IPROBLEM,'hessian_approximation',
     1'limited-memory')


      if (IERR.ne.0 ) goto 9990
C
C     Set an integer option
C
      IERR = IPADDINTOPTION(IPROBLEM, 'max_iter',150)
      if (IERR.ne.0 ) goto 9990

      !IERR = IPADDINTOPTION(IPROBLEM, 'print_level',0)
      !if (IERR.ne.0 ) goto 9990
C
C     Set a double precision option
C
      IERR = IPADDNUMOPTION(IPROBLEM, 'tol', 1.d-7)
      if (IERR.ne.0 ) goto 9990
C      IERR = IPADDSTROPTION(IPROBLEM,'derivative_test',
C     1'first-order')

      IERR = IPADDNUMOPTION(IPROBLEM, 'derivative_test_perturbation',
     1 5.d-5)

       IERR = IPADDNUMOPTION(IPROBLEM, 'derivative_test_tol',
     1 1.d-3)     
 
       IERR = IPADDNUMOPTION(IPROBLEM, 'point_perturbation_radius',
     1 1.d-4) 
        IERR = IPADDNUMOPTION(IPROBLEM,'expect_infeasible_problem_ctol',
     1 1.d-7)
C     Call optimization routine
      p=3d0
      !call fvm_simulate(X,d,DAT,IDAT,p)
      !call adjoint(X,d,GRAD,DAT,IDAT,p)
      ! X=0.4d0


      IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, DAT, IDAT)

      call stopClock()
      call toc()

      !CALL CPU_TIME(finish)
      !print *,finish,start
      !WRITE(*, '(" CPU time: ", F10.6, " sec.")') (finish - start)
      !call toc()


       read(*,*)
      open (unit = 2, file = "optimization_result")
      do I=1,N
      write (2,*) X(I)
      enddo
      close (2)
      deallocate (DAT,IDAT,DAT3,DAT4)
C     Output:
C
      if( IERR.eq.IP_SOLVE_SUCCEEDED ) then
         write(*,*)
         write(*,*) 'The solution was found.'
         write(*,*)
         write(*,*) 'The final value of the objective function is ',F
         write(*,*)
         write(*,*) 'The optimal values of X are:'
         write(*,*)
C         do i = 1, N
C           write(*,*) X(i)
C         enddo
         write(*,*)
         write(*,*) 'The error code is ',IERR
         write(*,*)
      endif
C
 9000 continue
C
C     Clean up
C
      call IPFREE(IPROBLEM)
      call pardiso_release(d)
      stop
C
 9990 continue
      write(*,*) 'Error setting an option'
      goto 9000



      end

