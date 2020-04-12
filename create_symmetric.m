clear all; clc;


load('optimization_result')


n_design_cells = size(optimization_result,1)*4;

figure
fortran=  reshape(optimization_result,[sqrt(n_design_cells)/2,sqrt(n_design_cells)/2])';

fortran = [(fortran),fliplr(fortran);flipud(fortran), fliplr(flipud(fortran))];

x=linspace(1,10,50);
y=linspace(1,10,50);
heatmap(fortran,'gridvisible','off')
title('optimal solution')


