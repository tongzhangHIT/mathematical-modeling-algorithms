function [out]=fitness(in)
x=in(:,1);
y=in(:,2);
out=-20*exp(-0.2*sqrt(0.5*(x.^2+y.^2)))-exp(0.5*(cos(2*pi*x)+cos(2*pi*y)))+22.71282;