figure
[x,y]=meshgrid(-3:0.5:3,-3:0.1:3);
z=peaks(x,y);
ribbon(y,z)
xlabel('x'),ylabel('y'),zlabel('z')

figure
ribbon(z)