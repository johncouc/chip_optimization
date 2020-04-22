%% this file contains a solver of the optimization problem using fmincon. 
%% the code works, but gets very slow above 60x60 variables. nonetheless, it gives good results,
%% solutions with very similar structure to that of the Fortran solver.
clear all; clc;
close all;

Dx = 0.01/2;
Dy = 0.01/2;
Dz = 0.001;
load('optimization_result');

n_design_cells = 20*20;
dx = Dx/sqrt(n_design_cells);
dy = Dy/sqrt(n_design_cells);

beq = 0.4*Dx*Dy;
%area_constraint = ones(1,400)*0.0005*0.0005;
area_constraint = ones(1,n_design_cells)*dx*dy;
Aeq = area_constraint;
v = 0.4*ones(1,n_design_cells)';
lb = zeros(n_design_cells,1);
ub = 1*ones(n_design_cells,1);

opts = optimoptions('fmincon','Algorithm','interior-point','GradObj','on','display','iter','MaxFunEvals',100);

%opts = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',1,'display','iter');

FUN = @(z) fvm_with_Gradient(z);
FUN2 = @(z) just_fvm(z);

% switch: material constraint equality/inequality
[v_min,FVAL,exitflag,output,lambda,grad,hessian]=fmincon(FUN,v,[],[],Aeq,beq,lb,ub,[],opts);

%[v_min,FVAL,exitflag,output,lambda,grad,hessian]=fmincon(FUN,v,Aeq,beq,[],[],lb,ub,[],opts);

figure
vopt=  reshape(v_min,[sqrt(n_design_cells),sqrt(n_design_cells)])';

vopt = [(vopt),fliplr(vopt);flipud(vopt), fliplr(flipud(vopt))];
heatmap(vopt,'gridvisible','off')
title('optimal solution')

[cost,T]= just_fvm(v_min)
Tsize=floor(sqrt(size(T,1)));

T_unfold= reshape(T,Tsize,Tsize)';

T_unfold = [(T_unfold),fliplr(T_unfold);flipud(T_unfold), fliplr(flipud(T_unfold))];

 figure
 heatmap(T_unfold,'gridvisible','off')
 title('temp distribution for optimal solution')
% plot_fvm_2(v_min)
 %save('optim_matlab','v_min')
% v=  reshape(v,[sqrt(n_design_cells),sqrt(n_design_cells)])';
% 
% v = [(v),fliplr(v);flipud(v), fliplr(flipud(v))];
% heatmap(v,'gridvisible','off')
% title('optimal solution')