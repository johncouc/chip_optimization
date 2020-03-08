! COMMANDS TO RUN:
! gfortran -c main2.f90 eigv_module.f90 fvm.f90
! gfortran gfortran -o main2 main2.o sparse.o fvm.o -llapack

! ./main2


PROGRAM main2

	USE sparse_system_solvers
	USE eigv_module
        USE fvm

	INTEGER, PARAMETER		:: d = 6
	real, dimension(5,5) :: v = reshape(&
(/0.4,0.4,0.4,0.4,0.4, &
0.4,0.4,0.4,0.4,0.4,&
0.4,0.4,0.4,0.4,0.4,&
0.4,0.4,0.4,0.4,0.4,&
0.4,0.4,0.4,0.4,0.4/), (/5,5/))


call fvm_with_gradient(v, d)

print *,"result", v

END PROGRAM
