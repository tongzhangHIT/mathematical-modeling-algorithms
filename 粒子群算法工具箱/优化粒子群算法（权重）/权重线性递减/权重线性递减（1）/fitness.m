function F = fitness(x)
F=100*(x(1)^2-x(2))^2+(1-x(1))^2;
%[xm,fv]= LinWPSO(@fitness,N,c1,c2,w1,w2,M,D)