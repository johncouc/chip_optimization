%% COMPARISON OF GRAD CALCULATED BY ADJOINT METHOD TO FMINCON'S APPROXIMATION

clear all; clc;
close all;
%load("A");
%load("T");
%load("f");


Dx = 0.01/2;
Dy = 0.01/2;
Dz = 0.001;


n_design_cells =10*10;
dx = Dx/sqrt(n_design_cells);
dy = Dy/sqrt(n_design_cells);

% Optimization of the objective function
beq = 0.4*Dx*Dy;
%area_constraint = ones(1,400)*0.0005*0.0005;
area_constraint = ones(1,n_design_cells)*dx*dy;
Aeq = area_constraint;
v = rand(n_design_cells,1);
%v=0.4*ones(n_design_cells,1);
lb = zeros(n_design_cells,1);
ub = 1000*ones(n_design_cells,1);

opts = optimoptions('fmincon','Algorithm','interior-point','GradObj','off','display','iter','MaxFunEvals',1,"FiniteDifferenceStepSize",10^(-5));
%opts = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',1,'display','iter');

FUN = @(z) fvm_with_Gradient(z);

FUN2 = @(z) just_fvm(z);



[v_min,FVAL,exitflag,output,lambda,grad,hessian]=fmincon(FUN,v,[],[],Aeq,beq,[],[],[],opts);


%AMATRIX = reshape(AMATRIX',121,121) 
[cost,grad_adj]=fvm_with_Gradient(v)
checker=abs((grad_adj'-grad)./grad)
 %solFin= reshape(grad_adj,[sqrt(n_design_cells),sqrt(n_design_cells)])';
 % heatmap(solFin);
  
  solFin2= reshape(checker,[sqrt(n_design_cells),sqrt(n_design_cells)])';
  figure
  heatmap(solFin2);
  title 'relative deviation'
  
  
% additional finite difference tests for isolated variables can be done
% here
step =zeros(n_design_cells,1);
eps=10^(-4);

step(10)=eps;
base =fvm_with_Gradient(v);
ble=v+step;
forward = fvm_with_Gradient(v+step);
deriv= (forward-base)/eps;
abs((grad_adj(10)-deriv)./deriv)
