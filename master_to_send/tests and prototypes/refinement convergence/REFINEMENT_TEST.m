%% EXAMPLE FILE SHOWING THE SOLUTION AFTER REFINEMENT

clear all; clc;
close all;
%load("A");
%load("T");
%load("f");


Dx = 0.01/2;
Dy = 0.01/2;
Dz = 0.001;


n_design_cells =20*20;
dx = Dx/sqrt(n_design_cells);
dy = Dy/sqrt(n_design_cells);

% Optimization of the objective function
beq = 0.4*Dx*Dy;
%area_constraint = ones(1,400)*0.0005*0.0005;
area_constraint = ones(1,n_design_cells)*dx*dy;
Aeq = area_constraint;
v = 18*rand(n_design_cells,1);
%v=0.4*ones(n_design_cells,1);
lb = zeros(n_design_cells,1);
ub = 1000*ones(n_design_cells,1);

opts = optimoptions('fmincon','Algorithm','interior-point','GradObj','off','display','iter','MaxFunEvals',1,"FiniteDifferenceStepSize",10^(-5));
%opts = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',1,'display','iter');

FUN = @(z) fvm_with_Gradient(z);

FUN2 = @(z) just_fvm(z);
    




%AMATRIX = reshape(AMATRIX',121,121) 
[cost,grad_adj]=fvm_with_Gradient(v)
 %solFin= reshape(grad_adj,[sqrt(n_design_cells),sqrt(n_design_cells)])';
 % heatmap(solFin);
    load('T_matlab')
    
 T_1 = T;
  solFin2= fliplr(reshape(T_1,floor(sqrt(n_design_cells)+1),floor(sqrt(n_design_cells)+1) ));
  figure
  heatmap(solFin2);
  grid off;
  title 'Temperature field'
  
  n_design_cells=4*n_design_cells;
  v2=refine(v);
  [cost,grad_adj]=fvm_with_Gradient(v2)
load('T_matlab')
    
 T_2 = T;
    solFin3= fliplr(reshape(T_2,floor(sqrt(n_design_cells)+1),floor(sqrt(n_design_cells)+1) ));
  figure
  heatmap(solFin3);
  title 'Temperature field - Refined'
  grid off;
  
  n_design_cells=4*n_design_cells;
  v3=refine(v2);
  [cost,grad_adj]=fvm_with_Gradient(v3)
load('T_matlab')
    
 T_3 = T;
    solFin4= fliplr(reshape(T_3,floor(sqrt(n_design_cells)+1),floor(sqrt(n_design_cells)+1) ));
  figure
  heatmap(solFin4);
  title 'Temperature field - Refined'
  grid off;
