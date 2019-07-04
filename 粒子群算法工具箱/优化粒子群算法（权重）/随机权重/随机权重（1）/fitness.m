function F =fitness(x)
F=-20*exp(-0.2*sqrt((x(1)^2+x(2)^2)/2))-exp((cos(2*pi*x(1))+cos(2*pi*x(2)))/2)+20+2.71928;
%[xm,fv] = RandWPSO(@fitness,N,c1,c2,mean_max,mean_min,sigma,M,D)
