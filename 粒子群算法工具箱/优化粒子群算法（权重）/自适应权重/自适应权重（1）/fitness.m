function F = fitness(x)

F=100*(x(1)^2-x(2))^2+(1-x(1))^2;
%[xm,fv] = SAPSO(@fitness,N,c1,c2,wmax,wmin,M,D)