C Copyright (C) 2002, 2010 Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Eclipse Public License.
C
C    $Id$
C
C =============================================================================
C
C     This is an example for the usage of IPOPT.
C     It implements problem 71 from the Hock-Schittkowski test suite:
C
C     min   x1*x4*(x1 + x2 + x3)  +  x3
C     s.t.  x1*x2*x3*x4                   >=  25
C           x1**2 + x2**2 + x3**2 + x4**2  =  40
C           1 <=  x1,x2,x3,x4  <= 5
C
C     Starting point:
C        x = (1, 5, 5, 1)
C
C     Optimal solution:
C        x = (1.00000000, 4.74299963, 3.82114998, 1.37940829)
C
C =============================================================================
C
C
C =============================================================================
C
C                            Main driver program
C
C =============================================================================
C
      program example

C
      use fvm
        
      implicit none
C
C     include the Ipopt return codes
C
      include 'IpReturnCodes.inc'
      
C     Size of the problem (number of variables and equality constraints)
C
      integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
      parameter  (N = 63*63, M = 1, NELE_JAC = N, NELE_HESS = 0)
      parameter  (IDX_STY = 1 )
       INTEGER, PARAMETER		:: d = 64
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
C     This could be used to pass double precision and integer arrays untouched
C     to the evaluation subroutines EVAL_*
C
      double precision, dimension (:,:), allocatable ::DAT
      double precision, dimension (:,:), allocatable:: IDAT
C
C     Place for storing the Ipopt Problem Handle
C
CC     for 32 bit platforms
C      integer IPROBLEM
C      integer IPCREATE
C     for 64 bit platforms:
      integer*8 IPROBLEM
      integer*8 IPCREATE
C
      integer IERR
      integer IPSOLVE, IPADDSTROPTION
      integer IPADDNUMOPTION, IPADDINTOPTION
      integer IPOPENOUTPUTFILE
C
      double precision F
      integer i
C
C     The following are the Fortran routines for computing the model
C     functions and their derivatives - their code can be found further
C     down in this file.
C
C      external  EV_G,  EV_JAC_G 
CC
CC     The next is an optional callback method.  It is called once per
CC     iteration.
CC
C      external ITER_CB

      
        
      real(8), dimension(1,N) :: v , GRAD
            !  call fvm_simulate(v, d,DAT,IDAT)
C     d**2,d**2+1
      allocate (DAT( d**2,d**2+1))
      allocate (IDAT (d+1,d+1))

      DAT=0
      IDAT=0 
      v=0.4d0 
 
  
     













C     Set initial point and bounds:
C
      X =  0.0
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
      IERR = IPADDINTOPTION(IPROBLEM, 'max_iter', maxiters)
      if (IERR.ne.0 ) goto 9990
C
C     Set a double precision option
C
      IERR = IPADDNUMOPTION(IPROBLEM, 'tol', 1.d-7)
      if (IERR.ne.0 ) goto 9990

CC
CC     Set a callback function to give you control once per iteration.
CC     You can use it if you want to generate some output, or to stop
CC     the optimization early.
CC
C      call IPSETCALLBACK(IPROBLEM, ITER_CB)

C
C     As a simple example, we pass the constants in the constraints to
C     the EVAL_C routine via the "private" DAT array.
C
C      DAT(1) = 0.d0
C      DAT(2) = 0.d0
C
C     Call optimization routine
      p=3d0
      IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, DAT, IDAT)
       read(*,*)
      open (unit = 2, file = "optimization_result")
      do I=1,N
      write (2,*) X(I)
      enddo
      close (2)
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
      stop
C
 9990 continue
      write(*,*) 'Error setting an option'
      goto 9000



      end


C =============================================================================
C
C                Computation of gradient of objective function
C
C =============================================================================
C
    
C
C =============================================================================
C
C                     Computation of equality constraints
C
C =============================================================================
C
     
C =============================================================================
C
C                Computation of Jacobian of equality constraints
C
C =============================================================================
C
C =============================================================================
C
C                Computation of Hessian of Lagrangian
C
C =============================================================================
C
C
C =============================================================================
C
C                   Callback method called once per iteration
C
C =============================================================================
C
      subroutine ITER_CB(ALG_MODE, ITER_COUNT,OBJVAL, INF_PR, INF_DU,
     1     MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT,
     2     DAT, ISTOP)
      implicit none
      integer ALG_MODE, ITER_COUNT, LS_TRIAL
      double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
      double precision ALPHA_DU, ALPHA_PR
      double precision DAT(*), IDAT(*)
      
      integer ISTOP
C
C     You can put some output here
C
      write(*,*) 'Testing callback function in iteration ', ITER_COUNT
C
C     And set ISTOP to 1 if you want Ipopt to stop now.  Below is just a
C     simple example.
C
      if (INF_PR.le.1D-04) ISTOP = 1

      return
      end
