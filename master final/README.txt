Year: 2019-2020
Project: heat 
Author 1: Robert Bogle
Author 2: Erikas Svazas
Author 3: Ioannis Koukoulis 
Needed libraries:
* Intel Math Kernel Library
* IPOpt 3.13.1. 
Instructions:

### DESCRIOTION ###
Program consists of 4 modules: 
* main.f90, 
* fvm.f90                      (contains subroutines for FVM simulation, adjoint, objective function and constraint evaluation),
* pardiso_sparse_solver_v3.f90 (contains subroutines and configuration of the sparse solver),
* timings.f90                  (used for timing tests).

### COMPILING/LINKING ###
to compile/link, you can look at the "compile" scripts, which is what we used on different operating systems.

### PROBLEM OPTIONS ###
In the main file you can set all options for the optimization. This includes options such as:
no of variables and no of iterations. One can set the initial value (or continuing where you left off by uncommenting lines).
There is a plethora of options that can be set for IPopt: options about line search, barrier strategy,
sparse solver options, statistics and ...
The settings that are already there, are the ones that worked best for our optimization
but there is definitely a room for exploration, especially in the commented options.
For more information, please consult the documentation: https://coin-or.github.io/Ipopt/OPTIONS.html

### OUTPUT ###
The output of the code is a file "optimization_result". This file contains 1/4 of the solution vector.
For a visualization of the final result, please run the Matlab script "create_symmetric.m".
The file "ipopt.out" contains the terminal output of the program.

### TESTING ###
* To run the derivative test, enable the the derivative-test option.
* One can change the SIMP strategy and have a varying penalty parameter through the optimization
  (default is p=3 constant).
* The sparse solver can be changed to iterative with preconditioning, by toggling MODE = 1 
  in the pardiso_sparse_solver_v3.f90 file.
The folder "tests" contains matlab files of prototyping and the code verification cases.


Can we use and modify your code for demonstrations concerning the course? Yes 
