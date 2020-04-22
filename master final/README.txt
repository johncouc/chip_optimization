Year: 2019-2020
Project: heat 
Author 1: Robert Bogle
Author 2: Erikas Svazas
Author 3: Ioannis Koukoulis 
Needed libraries:
* Intel Math Kernel Library
* IPOpt 3.13.1. 
Instructions:

Program consists of 4 modules : main.f90, fvm.f90 (contains subroutines for 
fvm simulation, adjoint, objective function and cnstrnt evaluation), pardiso_sparse_solver_v3.f90 (contains subroutines
and configuration of the sparse solver) timings.f90 (used for timing tests)

to compile/link, you can look at the "compile" scripts, which is what were used in the different environments

Using the main file you can set all options for the optimization. This includes basic options such as:
no of variables and no of iterations.
setting initial value (or continuing where you left off by uncommenting lines )

There's also a plethora of options that can be set for IPopt.
there are options about line search, barrier strategy, sparse solver options,statistics and so on.

The seetings that are already there, are the ones that seem to work best for our optimization
 but there's definitely room for exploration, especially in the commented options.
for more info, consult the documentation : https://coin-or.github.io/Ipopt/OPTIONS.htm



Output of the code is a file optimization_result. this contains 1/4 of the solution vector. you can run the matlab script
create_symmetric.m for a visualization of the final result.

The file ipopt.out is also created, containing the screen output of the program. 

* A Derivative test can be made by enabling the derivative-test option.
* You can change the SIMP strategy, so as to have a varying penalty parameter through the optimization (default is 
p =3 constant)
* The sparse solver can be changed to iterative with preconditioning, by toggling MODE = 1 
in the pardiso_sparse_solver_v3.f90 file.

the folder "tests" contains matlab files of prototyping and the code verification cases.


Can we use and modify your code for demonstrations concerning the course? Yes 
