Dx = 0.01
Dy = 0.01
Dz = 0.001

Ts= 293

n= 100 % number of subdomains ( dy should preferrably be a multiple of 0.001 ).

%% after discretization, it is decided how much metal to put in a particular cell.
v=ones(n,1)
k=zeros(sqrt(n),sqrt(n));
%% create discretized conductivity field (works for n= 100*a^2)
for i = 1:sqrt(n)
    for j= 1:sqrt(n)
    k(i,j)= 65*v((i-1)*sqrt(n)+j) + 0.2* ( 1 - v((i-1)*sqrt(n)+j)) ; % change this to SIMP
    end
end

%% constants 
a0=zeros(sqrt(n),sqrt(n));
a1=zeros(sqrt(n),sqrt(n));
a2=zeros(sqrt(n),sqrt(n));
a3=zeros(sqrt(n),sqrt(n));
a4=zeros(sqrt(n),sqrt(n));

for i = 2:sqrt(n)-1
    for j= 2:sqrt(n)-1
    ii = i - 1;
    jj = j - 1;
    a0(ii,jj)= k(i,j) + k(i-1,j) + k(i,j-1) + k(i-1,j-1);
    a1(ii,j)= 0.5* (k(i,j) + k(i,j-1) );
    a2(i,j)= 0.5* (k(i-1,j) + k(i,j) );
    a3(i,j)= 0.5* (k(i-1,j-1) + k(i-1,j) );
    a4(i,j)= 0.5* (k(i-1,j-1) + k(i,j-1) );
        end
end


