module fvm

	USE sparse_system_solvers
	USE eigv_module
	

IMPLICIT NONE
SAVE

CONTAINS

SUBROUTINE fvm_with_Gradient(v, d)
integer, intent(in)::d
	REAL				       :: l,b,a,i, j, ii, jj, m !number of design variables
	REAL                                   :: Cost, lengthx, lengthy, lengthz, Q, Ts, p
	REAL                                   ::  n, dx, dy, hx, hy, dkdk
	REAL                                   ::  bdir, dbdir, dSdir, Sdir
	REAL                                   :: ke, kw, kn, ks, dharm, c
	REAL				       :: dx_east, dx_west, dy_north, dy_south
	REAL, DIMENSION(3)                     :: dxi_dv, vec
	REAL, DIMENSION(5)                     :: xi
	REAL, DIMENSION((d-1)**2)                     :: G
	REAL, DIMENSION((d-1)**2), INTENT(inout)      :: v
	REAL, DIMENSION(d*d) :: f, PSI, theta_R_v, T,finv
	REAL, DIMENSION(d-1,d-1) 	   :: theta_kk_v, kk
	REAL, DIMENSION(d+1,d+1) :: theta_k_v, k
	REAL, DIMENSION(d*d,d*d) :: dxi_dk,  aa

	! assign variables
	
	lengthx = 0.01
	lengthy = 0.01
	lengthz = 0.001
	Q = 2/(lengthz*lengthx*lengthy)
	Ts = 293
	p = 3
	a = d
        b = a
        m = (a-1)**2
	n = a*b      !number of state grid cells. 
	dx = lengthx/(a-1) 
	dy = lengthy/(b-1)  !dy should preferrably divide 0.001
	theta_k_v = 0 
	dxi_dk = 0

	! after discretization, it is decided how much metal
    ! to put in a particular cell.
	


    kk = 1
	
    ! create discretized conductivity (design) field 
	
    do ii = 1,b-1
	do jj = 1,a-1
		kk(ii,jj) = 0.2+(65-0.2)*v(ii+(a-1)*(jj-1))**p  
		theta_kk_v(ii,jj) =  p*(65-0.2)*v(ii+(a-1)*(jj-1))**(p-1) 
	enddo
    enddo
	
	k = 0
    	k(2:a,2:b) = kk 
	k(1,2:a) = kk(1,:)
	k(a+1,2:a) = kk(a-1,:)     !fictitious design nodes 
	k(2:a,1) = kk(:,1)
	k(2:a,a+1) = kk(:,a-1)
	
   	theta_k_v(2:a,2:b) = theta_kk_v 
	theta_k_v(1,2:a) = theta_kk_v(1,:)
	theta_k_v(a+1,2:a) = theta_kk_v(a-1,:)     !fictitious design nodes 
	theta_k_v(2:a,1) = theta_kk_v(:,1)
	theta_k_v(2:a,a+1) = theta_kk_v(:,a-1)
	
	f = 0
	aa = 0

	do i = 1,n
		aa(i,i) = 1   
	enddo

print *, "A",aa,"a",a,"K",kk ,"k",k

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
    aa(l,l-1) =  xi(1) 
    dxi_dk(l,l-1) = -hy/dx_west 
  
end if

if (.not. (i==a)) then
    if (i==a-1 .or. i==1) then
    dx_east = 0.75*dx 
    else
    dx_east = dx 
    endif
    xi(2) =  -(ke*hy)/dx_east 
    aa(l,l+1) =  xi(2) 
    dxi_dk(l,l+1) =  -hy/dx_east 
end if

if (.not. (j==1)) then
    if (j==2 .or. j==b) then
    dy_south = 0.75*dy 
    else
    dy_south = dy 
    end if
    xi(3) =  - (ks*hx)/dy_south 
    aa(l,l-a) =  xi(3) 
    dxi_dk(l,l-a) = -hx/dy_south 
endif

if (.not. (j==b)) then
    if ((j==b-1) .or. (j==1)) then
    dy_north = 0.75*dy 
    else
    dy_north = dy 
    end if
 xi(4) = -(kn*hx)/dy_north 
 aa(l,l+a) = xi(4) 
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

aa(l,l) =  -(sum(xi))+(Sdir)
dxi_dk(l,l) = sum(dxi_dk(l,:))
f(l) = Q*hx*hy + bdir

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LINEAR SOLVE

!T = A\f
CALL solver_v1(aa, f, d*d, T)

!Gradient Calculations

!FIELD ADJOINT EQ:
!PSI = - T'*(A)^-1                      !good
!PSI=((A)\(-T))'                        !better
CALL solver_v1(aa,-T, d*d, PSI)         !best
!PSI = matmul(-T,aainv)                 !worst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LINEAR SOLVE

!print *, "A",aa,"T",T,"d",d ,"f", f ,"PSI",PSI

do l = 1, m
 
    theta_R_v = 0
    ii  =  1+ mod(l-1,a-1)
    jj =  ceiling((l)/(b-1)) 
    i =  1+ mod(l-1,a)
    j =  ceiling((l)/(b)) 
    m = ii+(a)*(jj-1)
    
	dxi_dv = 0

! 1st ELEMENT

dbdir = 0
dSdir = 0

if (jj==1) then
    dkdk = 1
else
    call harmonic_drv(k(ii+1,jj+1),k(ii+1,jj),dkdk)
endif

dxi_dv(1) = dxi_dk(m,m+1)*dkdk*theta_k_v(ii+1,jj+1)!east

if (ii==1) then
    dkdk = 1
    if (( ((jj-1)*dy >= 0.003)  .and. ((jj-1)*dy <= 0.007))) then
	call harmonic_drv(k(ii+1,jj+1),k(ii+1,jj), dharm)
    	call getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir)
    else
 	dbdir = 0
   	dSdir = 0
    end if
else
    call harmonic_drv(k(ii+1,jj+1),k(ii,jj+1),dkdk) 
endif

dxi_dv(2) = dxi_dk(m,m+a)*dkdk*theta_k_v(ii+1,jj+1)!north
dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2)) + dSdir
vec = (/T(m+1),T(m+a),T(m)/)
theta_R_v(m) =  dot_product(dxi_dv, vec) - dbdir

!! 2nd ELEMENT

dbdir = 0
dSdir = 0

if (jj==1) then
    dkdk = 1
else
    call harmonic_drv(k(ii+1,jj+1),k(ii+1,jj),dkdk)
end if

dxi_dv(1) = dxi_dk(m+1,m)*dkdk*theta_k_v(ii+1,jj+1)!west

if (ii==a-1) then
    dkdk = 1
    if (( ((jj-1)*dy >= 0.003)  .and. ((jj-1)*dy <= 0.007))) then
    call harmonic_drv(k(ii+1,jj+1),k(ii+1,jj), dharm)
    call getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir) 
    else
    dbdir =  0
    dSdir = 0
    end if
else
    call harmonic_drv(k(ii+1,jj+1),k(ii+2,jj+1),dkdk)
    dbdir =  0
    dSdir = 0 
end if
dxi_dv(2) = dxi_dk(m+1,m+a+1)*dkdk*theta_k_v(ii+1,jj+1)!north
dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2)) + dSdir
vec = (/T(m),T(m+1+a),T(m+1)/)
theta_R_v(m+1) =  dot_product(dxi_dv, vec) - dbdir
dbdir =  0
    dSdir = 0
!! 3rd ELEMENT
if (jj==b-1) then
    dkdk = 1
else
    call harmonic_drv(k(ii+1,jj+1),k(ii+1,jj+2),dkdk)
endif

dxi_dv(1) = dxi_dk(m+a,m+a+1)*dkdk*theta_k_v(ii+1,jj+1)!east
if (ii==1) then
    dkdk = 1
    if (( ((jj)*dy >= 0.003)  .and. ((jj)*dy <= 0.007)))  then
call harmonic_drv(k(ii+1,jj+1),k(ii+1,jj+2), dharm)
    call getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir) 	
     else
    dbdir =  0
    dSdir = 0
    endif
else
    call harmonic_drv(k(ii+1,jj+1),k(ii,jj+1),dkdk)
    dbdir = 0
    dSdir = 0
endif
dxi_dv(2) = dxi_dk(m+a,m)*dkdk*theta_k_v(ii+1,jj+1)!south
dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2))+dSdir
vec = (/T(m+a+1),T(m),T(m+a)/)
theta_R_v(m+a) =  dot_product(dxi_dv, vec) - dbdir

!! 4th ELEMENT

dbdir =  0
    dSdir = 0
if (jj==b-1) then
    dkdk = 1
else
    call harmonic_drv(k(ii+1,jj+1),k(ii+1,jj+2),dkdk)
endif

dxi_dv(1) = dxi_dk(m+a+1,m+a)*dkdk*theta_k_v(ii+1,jj+1)!west

if (ii==a-1) then
    dkdk = 1
    if (( ((jj)*dy >= 0.003)  .and. ((jj)*dy <= 0.007))) then
	call harmonic_drv(k(ii+1,jj+1),k(ii+1,jj+2), dharm)
    call getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir) 
    else
    dbdir =  0
    dSdir = 0
    endif
else
    call harmonic_drv(k(ii+1,jj+1),k(ii+2,jj+1),dkdk)
    dbdir = 0
    dSdir = 0
endif 
dxi_dv(2) = dxi_dk(m+a+1,m+1)*dkdk*theta_k_v(ii+1,jj+1)!south
dxi_dv(3) =  -(dxi_dv(1)+dxi_dv(2)) + dSdir
vec = (/T(m+a),T(m+1),T(m+a+1)/)
theta_R_v(m+a+1) =  dot_product(dxi_dv, vec) - dbdir
G(l) = dot_product(PSI,theta_R_v)
enddo
v = G
cost = 0.5* dot_product(T,T)
end subroutine

subroutine harmonic_drv(x,c,dharm) 
real, intent(out) :: dharm
real, intent(in) :: x,c
dharm= 2*(c**2/(x+c)**2)
end subroutine 



subroutine getdir(d,ii,jj,dy,dx,Ts,dharm,theta_k_v,dbdir,dSdir)
integer    , intent(in)             :: d
	REAL			:: ii, jj
	REAL                  :: dx, dy, Ts, dbdir, dSdir, dharm
	REAL, DIMENSION(d*d) :: PSI
	REAL, DIMENSION(d+1,d+1) :: theta_k_v

    dbdir =  4*Ts*dharm*theta_k_v(ii+1,jj+1)*dy/dx
    dSdir = 4*dharm*theta_k_v(ii+1,jj+1)*dy/dx  
end subroutine

END MODULE




