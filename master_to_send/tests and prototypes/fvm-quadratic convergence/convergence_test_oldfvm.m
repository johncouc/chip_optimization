clear
close all;
Dx = 0.01;
Dy = 0.01;
Dz = 0.001;
Q=2/(Dz*Dx*Dy);
Ts= 0;
for counter= 1:1
a=11*counter; %(elements/line)
b=11*counter; %(elements/col)
n= a*b; % number of state grid cells. 
N=(a-1)*(b-1);
dx=Dx/(a-1);
dy=Dy/(b-1); %( dy should preferrably divide 0.001 ).

%% after discretization, it is decided how much metal to put in a particular cell.
v=0*ones(N,1);
k=zeros(b-1,a-1);
%% create discretized conductivity (design) field 
for I = 1:b-1
    for J= 1:a-1
    k(I,J)= 1.65*v(J+(a-1)*(I-1)) + 0.2* ( 1 - v(J+(a-1)*(I-1)) ) ; % change this to SIMP
    
    end
end
k=[k(1,:);k;k(end,:)];
k =[k(:,1),k,k(:,end)];  %fictitious design nodes 

%% temperature field.corners: 1,a, n-a+1, n.  Dirichlet B.C:left at 1+(a*0.003/dy) - 1+(a*0.007/dy). 
%% Right at a + (a*0.003/dy) - a+(a*0.007/dy).
T = zeros(n,1);
f= Ts*ones(n,1);
A=eye(n,n);

for l=1:n
   i= 1+ mod(l-1,a);
   j= ceil((l)/b);
   hx=dx;
   hy=dy;
if (i==1) ;
    hx=dx/2;
    kw=0
    xcoord=hx/2;
else 
    kw= 0.5* (k(i,j)+k(i,j+1));
    xcoord=(i-1)*hx;
end
if (i==a)
    ke=0;
    hx=dx/2;
    xcoord=Dx-hx/2;
else
    ke = 0.5*(k(i+1,j)+k(i+1,j+1));
end
if (j==1)
    ks=0;
    hy=dy/2;
    ycoord=hy/2;
else
    ks = 0.5*(k(i,j) + k(i+1,j));
    ycoord=(j-1)*hy;
end
if (j==b)
    kn=0;
    hy=dy/2;
    ycoord=Dy-hy/2;
else
    kn = 0.5*(k(i,j+1) + k(i+1,j+1));
end
if ~(  ((i==1)|(i==a))   & ( ((j-1)*dy >= 0.000) & ((j-1)*dy <= 0.01))) 

if (l ~= 1)
A(l,l-1)=- (kw*hy)/dx;
end
if (l~= n)
A(l,l+1)= -(ke*hy)/dx;
end

if (j ~= b)
A(l,l+a)= - (kn*hx)/dy;
end

if (j~=1)
A(l,l-a)= - (ks*hx)/dy;
end    
A(l,l) = (ks+kn)*hx/dy +(ke+kw)*hy/dx;
f(l) = Q*hx*hy;
end

if (i==a) %test dirichlet cond.
f(l)= 300+30*cos(2*pi*ycoord/Dx);    
end
Tanal(l)=300/Dx*xcoord;
Tanal(l)=Tanal(l)+30/sinh(2*pi)*cos(2*pi*ycoord/Dx)*sinh(2*pi*xcoord/Dx);
Tanal(l)= Tanal(l)- Q/0.2*xcoord*(xcoord-Dx)/2;
    
end

T = A\f;
error(counter)=norm((T-Tanal'))/norm(Tanal);  % switch norms here
dh(counter)=dx;
end


loglog(dh,error)
hold on
loglog(dh,dh*15)
loglog(dh,(dh.^2)*35000)

legend('error curve','linear','quadratic')
hold off
solFin= reshape(T,[b,a])';
solAnalytic=reshape(Tanal,[b,a])'
ermap=reshape(T-Tanal',[b,a])';
heatmap(ermap);
heatmap(solFin);
title 'FVM solution'
figure
heatmap(solAnalytic)
title 'Analytic Solution'
%contourf((solFin))