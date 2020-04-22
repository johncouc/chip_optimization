%% Testing reduction on quarter of the chip. symmetric (random) v's produce symmetric T's.
%% the test is primitive (some preprocessing should be done to make sure 
%% the solutions match exactly), but shows that indeed, only 1/4th is needed.

clear all; clc;
close all;
;
Dx = 0.01;
Dy = 0.01;
Dz = 0.001;

n_design_cells =6*6;
dx = Dx/sqrt(n_design_cells);
dy = Dy/sqrt(n_design_cells);
d = sqrt(n_design_cells);
% Optimization of the objective function
beq = 0.4*Dx*Dy;
%area_constraint = ones(1,400)*0.0005*0.0005;
area_constraint = ones(1,n_design_cells)*dx*dy;
Aeq = area_constraint;
v=rand(n_design_cells,1);
%v=0.4*ones(n_design_cells,1);
lb = zeros(n_design_cells,1);
ub = 1000*ones(n_design_cells,1);

opts = optimoptions('fmincon','Algorithm','sqp','GradObj','off','display','iter','MaxFunEvals',1);
%opts = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',1,'display','iter');

FUN = @(z) fvm_with_Gradient_symm(z);

FUN2 = @(z) just_fvm(z);



%[v_min,FVAL,exitflag,output,lambda,grad,hessian]=fmincon(FUN,v,[],[],[],[],[],[],[],opts);
% for j = 1:10
% eps=(0.1)^(j)/size(v,1)
% for i= 1:size(v,1)
%     pert= zeros(size(v,1),1);
%     pert(i)=eps;%*sqrt(-1);
%     grad_fd(i,j)=(just_fvm(v+pert)-just_fvm(v))/(eps) ;
%     
% end
% end



[cost,grad_adj]=fvm_with_Gradient_symm(v)
 load('T_matlab')
T1=T;

T1 = reshape(T1,d+1,d+1);
T1=[T1,fliplr(T1);flipud(T1), fliplr(flipud(T1))]'

v_unfolded = reshape(v,d,d);
v_unfolded=[v_unfolded,fliplr(v_unfolded);flipud(v_unfolded), fliplr(flipud(v_unfolded))]
v_unfolded = reshape(v_unfolded,4*n_design_cells,1);

[cost,grad_adj]=fvm_with_Gradient(v_unfolded)
 load('T_matlab2')
 T2=T;
solFin= reshape(T2,(2*sqrt(n_design_cells)+1),(2*sqrt(n_design_cells)+1))';

%checker=abs((grad_adj'-grad)./grad)
 %solFin= reshape(,[sqrt(n_design_cells),sqrt(n_design_cells)])';
  heatmap(solFin);
  title('Solution on full chip')
  figure
  heatmap(T1);
  title('Solution on quarter of the chip')
 % solFin2= reshape(checker,[sqrt(n_design_cells),sqrt(n_design_cells)])';
 % figure
  %heatmap(solFin2);
  %title 'relative deviation'
  

% title 'FVM solution'