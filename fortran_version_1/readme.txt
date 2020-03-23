this is the 1st working version of the fortran solver for the cooling device problem.

It is based on the fvm.f90 and sparse_system_solvers.f90 files you guys made, combined with the example IPOPT file hs_071.f.

Changes in fvm.f90:

- fvm_with_gradient is now broken into two subroutines: FVM_SIMULATE(v,d,DAT,IDAT) and ADJOINT(v,d,GRAD,DAT,IDAT). 
The simulate part sets up the conductivities matrix IDAT,the stiffness matrix A and calculates the temperature field T.
 these last 2 are stored in the (d^2)x(d^2 + 1) matrix DAT. 

The adjoint part computes the theta_k_theta_v derivative from matrix IDAT, the adjoint temperature PSI from matrix DAT,
 sets up the rest of the derivatives and loops over all variables to get the gradient.


- There are 4 new subroutines, as used in the example code: EV_F_CHIP (calculates the actual cost based on the results of
fvm_simulate),  EV_GRAD_F_CHIP (fetches the gradient from the adjoint subroutine), EV_G (calculates the constraint, this
is the Aeq*v < beq part in the MATLAB code), and EV_JAC_G (calculates the jacobian of the constraint at each candidate
solution. normally this is a sparse matrix. in our case it is a non-sparse constant vector - seems to be working fine)

These 4 routines define the IPOPT problem in the main file(normally 1-2 additional routines for the hessian are required 
but we are using a BFGS method that approximates the hessian from the gradient, so these are not needed). 

Variables DAT and IDAT are used to communicate between these 4 subroutines at each solution, they were adapted from the example
problem code. fortran stupidly complains when trying to define their dimensions using the "d" integer, and this is why
i have to define them using these weird sqrt(floor functions. there is definitely a more elegant way to do it, but right now
it's working, i am tired of fortran's shit, so i don't care about fixing it.

other changes:

- the dxi_dk matrix has been changed from (d^2 x d^2) to d^2x5, since  there are at most 5 non-zero entries in each row.

-REAL was changed to REAL(8) ie made sure to be double precision.

-variables l,b,a,i, j, ii, jj, m were changed to their proper type, that is INTEGER (seriously what kind of sick MONSTER
declares i and j to be reals???)

-matrices kk, theta_kk were removed

-vector PSI was removed from getdir subroutine 


changes in sparse_system_solvers.f90:

-the code currently makes use of the classic LU decomposition subroutines. these have been changed from single to double precision
ie SGETRF was changed to DGETRF. 

-solver_v1 was modified: it does NOT copy the matrix anymore. our problem requires 2 solves, one for the temperature T and 
one forthe adjoint temperature PSI. these 2 solves have THE SAME left hand side, only the RHS is different. to do these  
with as little work as possible, a primitive fix was devised. Solver_v1 does the LU and stores the result in place of the original
matrix DAT(1:d,1:d). then solves for the temperature with DGETRS in DAT(1:d,d+1). The pivot vector IPIV is copied in a
 global variable IPIVV (CURRENTLY UP TO 99^2 VARS - CHANGE ITS SIZE IF YOU WANT TO SOLVE FOR GREATER # OF DESIGN VARIABLES).
The PSI vector is computed from the LU factorization that's already been performed, in subroutine solver_v11.

That way, instead of doing 2 copies and 2 decompositions of the A matrix, we make 0 copies and 1 decomposition.



 general :
-the program is compiled on my windows with the help of a (totally simplistic and unsophisticated) makefile called "batmake". this deletes all .o files
 and compiles everything from scratch, specifying all libraries used - by the way i haven't found a way to give the directory of the ipopt library files
so they have to be in the same directory as the source code. 

-the main hs_o71.f file is in FIXED form, while the other files are in free fortran form. this shouldn't be a problem, each program compiles fine on its own
and the linking is also fine. just know that when you try to copy-paste or comment out things you may get confusing errors.

-all arrays are defined on the heap using the -heap-arrays option, otherwise we can only solve for up to 400 vars or something.
-results are stored in files ipopt.out and solution vector in file 'scores'.
-even if you don't have ipopt you should still be able to compile the fvm and solver files. 


to do: 
-optimize everything! i will try and test the other solvers in the file. at some point it's best to convert the stiffness matrix in sparse format
and use the pardiso solver (less than 1% is nonzero entries for the # of design vars we're doing)

-I am curious to see if using parallel routines makes any difference - using sequential right now.

-definitely modify code to solve at 1/4th of the chip.

-get an even prettier solution figure. the one we have is 100x100 (takes about an hour), and looking very similar to the ones in the assignment pdf.
 but it  does look like an intermediate solution.

-use the simp method more methodically. currently the code uses p=3 from beginning to end, it's probably more useful to have it varying.

-changing the cost function from 0.5 TtT to 0.5 TtT/N^2 is best. this allows us to compare solutions of different fineness.

