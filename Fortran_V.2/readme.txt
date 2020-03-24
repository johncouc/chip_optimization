fortran v2

changes: 

-changed LU solver to cholesky decomposition routines DPOTRF DPOTRS (quite faster)

-took parameter p out of routines, it's sort of a global variable now. still haven't got a nice looking way to vary p within
the ipopt process. working with p=3 throughout the optimization works fine though. also added maxiters variable in fvm
(determines number of ipopt iterations)

-eliminated some unnecessary simulations happening at ev_f_chip

-made matrices data and idata ALLOCATABLE. they can become even larger now

-modified to solve for 1/4th of chip: Dx,Dy get halved, Q gets divided by 4, dirichlet boundary condition only applies for
x=0,  0.003<=y<=0.005. also constraint had to be slightly modified. 

-end result is stored in file 'optimization_result'. run the matlab code 'create_symmetric' to create a
 figure for the solution
