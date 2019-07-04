figure
rand('state',0);
x = rand(1,10);
y = rand(1,10);
TRI = delaunay(x,y);
triplot(TRI,x,y)
axis([0 1 0 1]);
hold on;
plot(x,y,'or');
hold off

numpts=192;
t=linspace(-pi,pi,numpts+1)';
r=0.1+5*sqrt(cos(6*t).^2+0.7^2);
x=r.*cos(t);
y=r.*sin(t);
dt=DelaunayTri(x,y);
tri=dt(:,:);
figure
triplot(tri,x,y)
axis equal