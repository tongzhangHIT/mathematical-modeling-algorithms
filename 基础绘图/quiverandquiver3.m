figure
[X,Y] = meshgrid(-2:.2:2); 
Z = X.*exp(-X.^2 - Y.^2);
[DX,DY] = gradient(Z,.2,.2);  %�ݶ�
contour(X,Y,Z)  %�ȸ���
hold on
quiver(X,Y,DX,DY)  %������
xlabel('x');
ylabel('y');
colormap hsv
hold off

figure
[X,Y] = meshgrid(-2:0.25:2,-1:0.2:1);
Z = X.* exp(-X.^2 - Y.^2);
[U,V,W] = surfnorm(X,Y,Z);  
quiver3(X,Y,Z,U,V,W,0.5);  %����
hold on
surf(X,Y,Z);
colormap hsv
view(-35,45)
axis ([-2 2 -1 1 -.6 .6])
hold off

vx = 2;     % x�����ٶ�
vy = 3;     % y�����ٶ�
vz = 10;    % z�����ٶ�
a = -32;    % z������ٶ�
 t = 0:.1:1; % ʱ��
x = vx*t;   % x����λ��
y = vy*t;   % y����λ��
z = vz*t + 1/2*a*t.^2;  % z����λ��
 u = gradient(x);        % x�����ݶ�
v = gradient(y);        % y�����ݶ�
w = gradient(z);        % z�����ݶ�
 scale = 0;
quiver3(x, y, z, u, v, w, scale)
view([70 18])

[x,y,z]=meshgrid(-0.8:0.2:0.8,-0.8:0.2:0.8,-0.8:0.2:0.8);
u=sin(pi*x).*cos(pi*y).*cos(pi*z);
v=-cos(pi*x).*sin(pi*y).*cos(pi*z);
w=sqrt(2/3)*cos(pi*x).*cos(pi*y).*sin(pi*z);
figure
quiver3(x,y,z,u,v,w)