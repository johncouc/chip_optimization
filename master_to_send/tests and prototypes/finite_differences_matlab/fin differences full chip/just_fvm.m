%% COMPUTES DIRICHLET BOUNDARY INDIRECTLY -GIVES QUADRATIC CONVERGENCE
function [Cost,g ]= just_fvm(v)
Dx = 0.01;
Dy = 0.01;
Dz = 0.001;
Q=2/(Dz*Dx*Dy);
Ts= 293;
p=3;
a=sqrt(size(v,1))+1; %(elements/line)
b=sqrt(size(v,1))+1; %(elements/col)
n= a*b; % number of state grid cells. 
N=(a-1)*(b-1);
dx=Dx/(a-1);
dy=Dy/(b-1); %( dy should preferrably divide 0.001 ).

%% after discretization, it is decided how much metal to put in a particular cell.

k=ones(b-1,a-1);
%% create discretized conductivity (design) field 
for I = 1:b-1
    for J= 1:a-1
    k(I,J)= 0.2 + (65-0.2)*v(J+(a-1)*(I-1))^p ;
    
    end
end
k=[k(1,:);k;k(end,:)];
k =[k(:,1),k,k(:,end)];  %fictitious design nodes 

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
%% calculate x,y coordinates for analytical solution:
if (i==1) ;
    xcoord=hx/2;
else
    xcoord=(i-1)*hx;
end
if (i==a)
    xcoord=Dx-hx/2;
end    
if (j==1)
    ycoord=hy/2;
else
    ycoord=(j-1)*hy;
end
if (j==b)
   ycoord=Dy-hy/2;
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
  
end

if ~(i==a) ;
    if (i==a-1|i==1)
    dx_east=0.75*dx;
    else
    dx_east=dx;
    end
    xi(2)= -(ke*hy)/dx_east;
    A(l,l+1)= xi(2);
end

if ~(j==1)
    if (j==2|j==b)
    dy_south=0.75*dy;
    else
    dy_south=dy;
    end
    xi(3)= - (ks*hx)/dy_south;
    A(l,l-a)= xi(3);
end

if ~(j==b)
    if (j==b-1|j==1)
    dy_north=0.75*dy;
    else
    dy_north=dy;
    end
 xi(4)= - (kn*hx)/dy_north;
 A(l,l+a)= xi(4);
end
%% define source terms for Dirichlet Boundary
if (( ((j-1)*dy >= 0.003) & ((j-1)*dy <= 0.007))) 
  
if (i == 1)
Ts=293 ;
bdir=2* Ts*kw*hy/hx;
Sdir=2*kw*hy/hx ;
end
if (i== a)
Ts=293     ; 
bdir= 2*Ts*ke*hy/hx;
Sdir=2*ke*hy/hx ;
end

end
A(l,l)= -(sum(xi))+(Sdir);
f(l) = Q*hx*hy + bdir;

end

T = A\f;
Cost=0.5*T'*T%/(norm(Ts*ones(n,1)));
end


% loglog(dh,error*9.5)
% hold on
% loglog(dh,dh*500)
% loglog(dh,(dh.^2)*350000)
% 
% legend('error curve','linear','quadratic')
% 
% hold off
% ermap=reshape(T-Tanal,[b,a])';
% heatmap(ermap);
% 
% solFin= reshape(T,[b,a])';
% solAnalytic=reshape(Tanal,[b,a])'
% 
% heatmap(solFin);
% title 'FVM solution'
% figure
% heatmap(solAnalytic)
% title 'Analytic Solution'


%contourf((solFin))