module FVM

IMPLICIT NONE
SAVE

CONTAINS

subroutine fvm_with_Gradient(v, m)

	INTEGER								   :: m, l !number of design variables
	REAL                                   :: Cost, Dx, Dy, Dz, Q, Ts, p
	REAL                                   :: a, b, n, dx, dy, hx, hy, dkdk
	REAL                                   :: i, j, I, J, l, bdir, dbdir, dSdir, Sdir
	REAL                                   :: ke, kw, kn, ks, dharm
	REAL								   :: dx_east, dx_west, dy_north, dy_south
	REAL, DIMENSION(3)                     :: dxi_dv, vec
	REAL, DIMENSION(5)                     :: xi
	REAL, DIMENSION(m)                     :: G
	REAL, DIMENSION(m), INTENT(inout)      :: v
	REAL, DIMENSION(sqrt(m)+1)*(sqrt(m)+1) :: f, PSI, theta_R_v, T
	REAL, DIMENSION(sqrt(m),sqrt(m)) 	   :: theta_K_v, K
	REAL, DIMENSION(sqrt(m)+2),sqrt(m)+2) :: theta_k_v, k
	REAL, DIMENSION(sqrt(m)+1)*(sqrt(m)+1,sqrt(m)+1)*(sqrt(m)+1) :: dxi_dk, f, A

	! assign variables
	
	Dx = 0.01
	Dy = 0.01
	Dz = 0.001
	Q = 2/(Dz*Dx*Dy)
	Ts = 293
	p = 3
	a = sqrt(m)+1  !(elements/line)
	b = sqrt(m)+1  !(elements/col)
	n = a*b      !number of state grid cells. 
	dx = Dx/(a-1) 
	dy = Dy/(b-1)  !dy should preferrably divide 0.001
	theta_k_v = 0 
	dxi_dk = 0

	! after discretization, it is decided how much metal
    ! to put in a particular cell.
	
    K = 1
	
    ! create discretized conductivity (design) field 
	
    do I = 1,b-1
		do J = 1,a-1
			K(I,J) = 0.2+(65-0.2)*v(I+(a-1)*(J-1))**p  
			theta_k_v(I,J) =  p*(65-0.2)*v(I+(a-1)*(J-1))**(p-1) 
		enddo
    enddo
	
	k = 0
    k(2:a,2:b) = K 
	k(1,2:a) = K(1,:)
	k(a+1,2:a) = K(a-1,:)     !fictitious design nodes 
	k(2:a,1) = K(:,1)
	k(2:a,a+1) = K(:,a-1)
	
    theta_k_v(2:a,2:b) = theta_K_v 
	theta_k_v(1,2:a) = theta_K_v(1,:)
	theta_k_v(a+1,2:a) = theta_K_v(a-1,:)     !fictitious design nodes 
	theta_k_v(2:a,1) = theta_K_v(:,1)
	theta_k_v(2:a,a+1) = theta_K_v(:,a-1)
	
	f = 0
	A = 0

	do i = 1,n
		A(i,i) = 1   
	enddo

do l = 1,n

   i= 1+ mod(l-1,a) 
   j =  ceiling(l/b) 
    
! check cell dimensions 

if ( i==a .or. i==1)  then
   hx = dx/2 
else 
   hx = dx 
endif

if (j==1 .or. j==b) then
   hy = dy/2 
else
   hy = dy 
endif

! initialize xi & interpolate for conductivities of each side

    kw =  2* (k(i,j)**(-1)+k(i,j+1)**(-1))**(-1) 
    ke =  2*(1/k(i+1,j)+1/k(i+1,j+1))**(-1) 
    kn = 2*(1/k(i,j+1) + 1/k(i+1,j+1))**(-1) 
    ks = 2*(1/k(i,j) + 1/k(i+1,j))**(-1)   
    xi = 0
    bdir = 0 
    Sdir = 0 
	
! start building stiffness matrix       

if (.not. (i==1)) then
    if (i==a .or. i==2) then
    dx_west = 0.75*dx 
    else
    dx_west = dx 
    end if
    xi(1) =  - (kw*hy)/dx_west 
    A(l,l-1) =  xi(1) 
    dxi_dk(l,l-1) = -hy/dx_west 
  
end if

if (.not. (i==a)) then
    if (i==a-1 .or. i==1) then
    dx_east = 0.75*dx 
    else
    dx_east = dx 
    end
    xi(2) =  -(ke*hy)/dx_east 
    A(l,l+1) =  xi(2) 
    dxi_dk(l,l+1) =  -hy/dx_east 
end if

if (.not. (j==1)) then
    if (j==2 .or. j==b) then
    dy_south = 0.75*dy 
    else
    dy_south = dy 
    end if
    xi(3) =  - (ks*hx)/dy_south 
    A(l,l-a) =  xi(3) 
    dxi_dk(l,l-a) = -hx/dy_south 
end

if (.not. (j==b)) then
    if ((j==b-1) .or. (j==1)) then
    dy_north = 0.75*dy 
    else
    dy_north = dy 
    end if
 xi(4) = -(kn*hx)/dy_north 
 A(l,l+a) = xi(4) 
 dxi_dk(l,l+a) = -hx/dy_north
 
end if

! define source terms for Dirichlet Boundary

if (( ((j-1)*dy >= 0.003) .and. ((j-1)*dy <= 0.007))) then
if (i == 1) then
Ts = 293 
bdir = 2* Ts*kw*hy/hx
Sdir = 2*kw*hy/hx 
end if
if (i== a) then
Ts = 293      
bdir =  2*Ts*ke*hy/hx
Sdir = 2*ke*hy/hx 
end if
endif

A(l,l) =  -(sum(xi))+(Sdir)
dxi_dk(l,l) = sum(dxi_dk(l,:))
f(l) = Q*hx*hy + bdir

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LINEAR SOLVE

!T = A\f
CALL solver_v1(A, f, (sqrt(m)+1), T)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LINEAR SOLVE

! Gradient Calculations

! FIELD ADJOINT EQ:

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LINEAR SOLVE

!PSI = - T'*(A)^-1                      !good
!PSI=((A)\(-T))'                        !better
CALL solver_v1(A,T, (sqrt(m)+1), PSI)   !best

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LINEAR SOLVE

!for each vi calculate theta R/ theta vi (1xn matrix): find the state cells(equations)
!which are affected by vi.
!if I,J are the coordinates of vi, affected cells are going to be (in state
!cell coordinates): (I,J) , (I+1,J), (I,J+1), (I+1,J+1). different xis are affected
!for each of those. (xi2,x4,x5), (xi1,xi4,xi5), (xi2,xi3,xi5) and
!(xi1,xi3,xi5) respectively. use chain rule, dot product with corresponding
!temps gives the values @ theta_R_theta_vi(no. of cell) 
!

 do l = 1, m
 
    theta_R_v = 0
    I  =  1+ mod(l-1,a-1)
    J =  ceiling((l)/(b-1)) 
    i =  1+ mod(l-1,a)
    j =  ceiling((l)/(b)) 
    m = I+(a)*(J-1)
    
dxi_dv = 0

! 1st ELEMENT

dbdir = 0
dSdir = 0

if (J==1) then
    dkdk = 1
else
    dkdk = harmonic_drv(k(I+1,J+1),k(I+1,J))
endif

dxi_dv(1) = dxi_dk(m,m+1)*dkdk*theta_k_v(I+1,J+1)!east

if (I==1) then
    dkdk = 1
    if (( ((J-1)*dy >= 0.003)  .and. ((J-1)*dy <= 0.007))) then
	
call harmonic_drv(k(I+1,J+1),k(I+1,J), dharm)
    dbdir =  4*Ts*dharm*theta_k_v(I+1,J+1)*dy/dx
    dSdir = 4*dharm*theta_k_v(I+1,J+1)*dy/dx  
    else
    dbdir = 0
    dSdir = 0
    end if
else
    dkdk = harmonic_drv(k(I+1,J+1),k(I,J+1)) 
endif

dxi_dv(2) = dxi_dk(m,m+a)*dkdk*theta_k_v(I+1,J+1)!north
dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2)) + dSdir
vec = (/T(m+1),T(m+a),T(m)/)
theta_R_v(m) =  dot_product(dxi_dv, vec) - dbdir

!! 2nd ELEMENT

dbdir = 0
dSdir = 0

if (J==1) then
    dkdk = 1
else
    dkdk = harmonic_drv(k(I+1,J+1),k(I+1,J))
end if

dxi_dv(1) = dxi_dk(m+1,m)*dkdk*theta_k_v(I+1,J+1)!west

if (I==a-1) then
    dkdk = 1
    
    if (( ((J-1)*dy >= 0.003)  .and. ((J-1)*dy <= 0.007))) then
	
call harmonic_drv(k(I+1,J+1),k(I+1,J), dharm)
    dbdir =  4*Ts*dharm*theta_k_v(I+1,J)*dy/dx
    dSdir = 4*dharm*theta_k_v(I+1,J)*dy/dx 
     else
    dbdir =  0
    dSdir = 0
    end if
else
    dkdk = harmonic_drv(k(I+1,J+1),k(I+2,J+1))
    dbdir =  0
    dSdir = 0 
end if
dxi_dv(2) = dxi_dk(m+1,m+a+1)*dkdk*theta_k_v(I+1,J+1)!north
dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2)) + dSdir
vec = (/T(m),T(m+1+a),T(m+1)/)
theta_R_v(m+1) =  dot_product(dxi_dv, vec) - dbdir
dbdir =  0
    dSdir = 0
!! 3rd ELEMENT
if (J==b-1) then
    dkdk = 1
else
    dkdk = harmonic_drv(k(I+1,J+1),k(I+1,J+2))
endif

dxi_dv(1) = dxi_dk(m+a,m+a+1)*dkdk*theta_k_v(I+1,J+1)!east
if (I==1) then
    dkdk = 1
    if (( ((J)*dy >= 0.003)  .and. ((J)*dy <= 0.007)))  then
call harmonic_drv(k(I+1,J+1),k(I+1,J+2), dharm)
    dbdir =  4*Ts*dharm*theta_k_v(I+1,J+1)*dy/dx
    dSdir = 4*dharm*theta_k_v(I+1,J+1)*dy/dx  	
     else
    dbdir =  0
    dSdir = 0
    endif
else
    call harmonic_drv(k(I+1,J+1),k(I,J+1),dkdk)
    dbdir = 0
    dSdir = 0
endif
dxi_dv(2) = dxi_dk(m+a,m)*dkdk*theta_k_v(I+1,J+1)!south
dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2))+dSdir
vec = (/T(m+a+1),T(m),T(m+a)/)
theta_R_v(m+a) =  dot_product(dxi_dv, vec) - dbdir

!! 4th ELEMENT
dbdir =  0
    dSdir = 0
if (J==b-1) then
    dkdk = 1
else
    dkdk = harmonic_drv(k(I+1,J+1),k(I+1,J+2))
endif

dxi_dv(1) = dxi_dk(m+a+1,m+a)*dkdk*theta_k_v(I+1,J+1)!west

if (I==a-1) then
    dkdk = 1
    
    if (( ((J)*dy >= 0.003)  .and. ((J)*dy <= 0.007))) then
	call harmonic_drv(k(I+1,J+1),k(I+1,J+2), dharm)
    dbdir =  4*Ts*dharm*theta_k_v(I+1,J+1)*dy/dx
    dSdir = 4*dharm*theta_k_v(I+1,J+1)*dy/dx     
    else
    dbdir =  0
    dSdir = 0
    endif
else
    call harmonic_drv(k(I+1,J+1),k(I+2,J+1),dkdk)
    dbdir = 0
    dSdir = 0
endif 

dxi_dv(2) = dxi_dk(m+a+1,m+1)*dkdk*theta_k_v(I+1,J+1)!south
dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2)) + dSdir
vec = (/T(m+a),T(m+1),T(m+a+1)/)
theta_R_v(m+a+1) =  dot_product(dxi_dv, vec) - dbdir

G(l) = dot_product(PSI,theta_R_v)

 enddo

Cost = 0.5* dot_product(T,T)

end subroutine

subroutine harmonic_drv(x,c,dharm) 
real, intent(out) :: dharm
real, intent(in) :: x,c
dharm= 2*(c^2/(x+c)^2);
end subroutine 
