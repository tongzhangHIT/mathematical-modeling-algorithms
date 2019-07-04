function F=fitness(x)
F=0.5+(sin(x(1)^2+x(2)^2)^2-0.5)/(1.0+0.001*(x(1)^2+x(2)^2))^2;
%[xm,fv] = SAPSO(@fitness,N,c1,c2,wmax,wmin,M,D)
