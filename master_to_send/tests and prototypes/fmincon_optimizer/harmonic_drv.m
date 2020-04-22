%% derivative of function f = 2/ (1/x + c)
function [dharm]= harmonic_drv(x,c) 
dharm= 2*(c^2/(x+c)^2);
end

