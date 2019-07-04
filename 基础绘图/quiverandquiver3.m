figure
[X,Y] = meshgrid(-2:.2:2); 
Z = X.*exp(-X.^2 - Y.^2);
[DX,DY] = gradient(Z,.2,.2);  %梯度
contour(X,Y,Z)  %等高线
hold on
quiver(X,Y,DX,DY)  %向量场
xlabel('x');
ylabel('y');
colormap hsv
hold off

figure
[X,Y] = meshgrid(-2:0.25:2,-1:0.2:1);
Z = X.* exp(-X.^2 - Y.^2);
[U,V,W] = surfnorm(X,Y,Z);  
quiver3(X,Y,Z,U,V,W,0.5);  %生成
hold on
surf(X,Y,Z);
colormap hsv
view(-35,45)
axis ([-2 2 -1 1 -.6 .6])
hold off

vx = 2;     % x方向速度
vy = 3;     % y方向速度
vz = 10;    % z方向速度
a = -32;    % z方向加速度
 t = 0:.1:1; % 时间
x = vx*t;   % x方向位移
y = vy*t;   % y方向位移
z = vz*t + 1/2*a*t.^2;  % z方向位移
 u = gradient(x);        % x方向梯度
v = gradient(y);        % y方向梯度
w = gradient(z);        % z方向梯度
 scale = 0;
quiver3(x, y, z, u, v, w, scale)
view([70 18])

[x,y,z]=meshgrid(-0.8:0.2:0.8,-0.8:0.2:0.8,-0.8:0.2:0.8);
u=sin(pi*x).*cos(pi*y).*cos(pi*z);
v=-cos(pi*x).*sin(pi*y).*cos(pi*z);
w=sqrt(2/3)*cos(pi*x).*cos(pi*y).*sin(pi*z);
figure
quiver3(x,y,z,u,v,w)