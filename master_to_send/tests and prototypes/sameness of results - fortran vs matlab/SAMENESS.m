%% EXAMPLE FILE COMPARING THE 2 VERSIONS' TEMP AND GRAD FOR A RANDOM DESIGN VECTOR.

%% TEMPERATURE & GRADIENT CALCULATED BY FORTRAN  : "temp" and "symmetric_grad"
clear all; clc;
close all;
load('optimization_result');
load("temp");
load("symmetric_grad");
Dx = 0.01/2;
Dy = 0.01/2;
Dz = 0.001;


n_design_cells =20*20;
dx = Dx/sqrt(n_design_cells);
dy = Dy/sqrt(n_design_cells);


v = optimization_result;
%

[cost,grad_adj]=fvm_with_Gradient(v)

  load('T_matlab')

error_T = max(abs(temp-T)) 
rel_T = max(abs(temp-T)./T)
error_G = max(abs(grad_adj'-symmetric_grad))
rel_G = max(abs(grad_adj- symmetric_grad')./grad_adj)
%% MAX ABS DEVIATION FOR TEMP VECTOR: 2e-11. , RELATIVE DEVIATION ~ e-14
%% MAX ABS DEVIATION FOR G VECTOR: 5e-09., RELATIVE DEVIATION ~e-16
%% in addition, the  forward relative error for the
%% linear system was computed (order of machine precision).