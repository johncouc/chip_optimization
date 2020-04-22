%% file used to test finite differences for full-chip adjoint

clear all; clc;
close all;
%load("A");
%load("T");
%load("f");
Dx = 0.01;
Dy = 0.01;
Dz = 0.001;
%load('v_found');
%v_found=v_min;
%opts = optimset('fmincon');
%opts.Algorithm = 'sqp';
%format long g
%[X,fval] = fmincon(@fvm2,[0.4],[],[],[],[],[],[],[],opts)

%%% ----- EXAMPLE -----
%numbPMF=3; %number  probability mass function (PMF) values
%p0=rand(1,numbPMF);p0=p0/sum(p0); %create a random PMF

%funOpt=@(p)(sum((p-p0).^2+(p-p0).^3)); %example of a nonlinear function
%pmf0=zeros(1,numbPMF); %initial guess of PMF
%lb = zeros(1,numbPMF); %lower bounds for probabilities
%ub = ones(1,numbPMF); %upper bounds for probabilities

%Aeq = ones(1,numbPMF); 
%beq=1;
%[pmfMin,F]=fmincon(funOpt,pmf0,[],[],Aeq,beq,lb,ub,[]);

n_design_cells =20*20;
dx = Dx/sqrt(n_design_cells);
dy = Dy/sqrt(n_design_cells);

% Optimization of the objective function
beq = 0.4*Dx*Dy;
%area_constraint = ones(1,400)*0.0005*0.0005;
area_constraint = ones(1,n_design_cells)*dx*dy;
Aeq = area_constraint;
v =rand(n_design_cells,1);
v=0.4*ones(n_design_cells,1);
lb = zeros(n_design_cells,1);
ub = 1000*ones(n_design_cells,1);

opts = optimoptions('fmincon','Algorithm','sqp','GradObj','off','display','iter','MaxFunEvals',1);
%opts = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',1,'display','iter');

FUN = @(z) fvm_with_Gradient(z);

FUN2 = @(z) just_fvm(z);



[v_min,FVAL,exitflag,output,lambda,grad,hessian]=fmincon(FUN,v,[],[],[],[],[],[],[],opts);
% for j = 1:10
% eps=(0.1)^(j)/size(v,1)
% for i= 1:size(v,1)
%     pert= zeros(size(v,1),1);
%     pert(i)=eps;%*sqrt(-1);
%     grad_fd(i,j)=(just_fvm(v+pert)-just_fvm(v))/(eps) ;
%     
% end
% end
load('AMATRIX')

AMATRIX = reshape(AMATRIX',121,121) 
[cost,grad_adj]=fvm_with_Gradi  ent(v)
checker=abs((grad_adj'-grad)./grad)
 %solFin= reshape(,[sqrt(n_design_cells),sqrt(n_design_cells)])';
  %heatmap(solFin);
  
  solFin2= reshape(checker,[sqrt(n_design_cells),sqrt(n_design_cells)])';
  figure
  heatmap(solFin2);
  title 'relative deviation'
  

% title 'FVM solution'