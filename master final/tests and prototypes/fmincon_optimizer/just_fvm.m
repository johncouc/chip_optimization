%% 1st Attempt at adjoint
function [Cost,T]= just_fvm(v)
Dx = 0.01/2;
Dy = 0.01/2;
Dz = 0.001;
Q=0.5/(Dz*Dx*Dy);
Ts= 293;
p=1;
a=sqrt(size(v,1))+1; %(elements/line)
b=sqrt(size(v,1))+1; %(elements/col)
n= a*b; % number of state grid cells. 
N=(a-1)*(b-1); % number of design variables
dx=Dx/(a-1);
dy=Dy/(b-1); %( dy should preferrably divide 0.001 ).

%theta_f_t=zeros(n,n);
%theta_A_v=zeros(n,N);
%theta_f_v=zeros(n,N);
theta_k_v=zeros(b-1,a-1);
dxi_dk=zeros(n,n);
%Psi= zeros(N,1);

%% after discretization, it is decided how much metal to put in a particular cell.

k=ones(b-1,a-1);
%% create discretized conductivity (design) field 
for I = 1:b-1
    for J= 1:a-1
    k(I,J)= 0.2 + (65-0.2)*v(I+(a-1)*(J-1))^p ;
    theta_k_v(I,J)= p*(65-0.2)*v(I+(a-1)*(J-1))^(p-1);
    end
end
k=[k(1,:);k;k(end,:)];
k =[k(:,1),k,k(:,end)];  %fictitious design nodes 
theta_k_v=[theta_k_v(1,:);theta_k_v;theta_k_v(end,:)];
theta_k_v =[theta_k_v(:,1),theta_k_v,theta_k_v(:,end)]; 


T = zeros(n,1);

f= 0*ones(n,1);
A=eye(n,n);

for l=1:n
   i= 1+ mod(l-1,a);
   j= ceil((l)/b);
    
%% check cell dimensions %%   
if (i==a|i==1)  
   hx=dx/2;
else   
   hx=dx;
end
if (j==1|j==b)
   hy=dy/2;
else
   hy=dy;
end   

%% initialize xi & interpolate for conductivities of each side
       %kw= 0.5* (k(i,j)+k(i,j+1));
       kw= 2* (k(i,j)^(-1)+k(i,j+1)^(-1))^(-1);
       
       
       %ke= 0.5*(k(i+1,j)+k(i+1,j+1));
       ke= 2*(1/k(i+1,j)+1/k(i+1,j+1))^(-1);
       %kn = 0.5*(k(i,j+1) + k(i+1,j+1));
       kn = 2*(1/k(i,j+1) + 1/k(i+1,j+1))^(-1);
       ks = 2*(1/k(i,j) + 1/k(i+1,j))^(-1);  
       %ks = 0.5*(k(i,j) + k(i+1,j));  
       xi= 0* ones(5,1);
       bdir=0;
       Sdir=0;
%% start building stiffness matrix       
if ~(i==1);
    if (i==a|i==2)
    dx_west=0.75*dx;
    else
    dx_west=dx;
    end
    xi(1)= - (kw*hy)/dx_west;
    A(l,l-1)= xi(1);
    dxi_dk(l,l-1)=-hy/dx_west;
  
end

if ~(i==a) ;
    if (i==a-1|i==1)
    dx_east=0.75*dx;
    else
    dx_east=dx;
    end
    xi(2)= -(ke*hy)/dx_east;
    A(l,l+1)= xi(2);
    dxi_dk(l,l+1)= -hy/dx_east;
end

if ~(j==1)
    if (j==2|j==b)
    dy_south=0.75*dy;
    else
    dy_south=dy;
    end
    xi(3)= - (ks*hx)/dy_south;
    A(l,l-a)= xi(3);
    dxi_dk(l,l-a)=-hx/dy_south;
end

if ~(j==b)
    if (j==b-1|j==1)
    dy_north=0.75*dy;
    else
    dy_north=dy;
    end
 xi(4)= - (kn*hx)/dy_north;
 A(l,l+a)= xi(4);
 dxi_dk(l,l+a)=-hx/dy_north;
end
%% define source terms for Dirichlet Boundary
if (( ((j-1)*dy >= 0.003) & ((j-1)*dy <= 0.005))) 
  
if (i == 1)
Ts=293 ;
bdir=2* Ts*kw*hy/hx;
Sdir=2*kw*hy/hx ;

end


end
A(l,l)= -(sum(xi))+(Sdir);
dxi_dk(l,l)= sum(dxi_dk(l,:) );
f(l) = Q*hx*hy + bdir;


end


T = A\f;
Cost=T'*T/(2*N)
end
