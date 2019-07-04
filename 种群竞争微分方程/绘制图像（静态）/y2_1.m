function y
%���Ƶ�r1=1;r2=1;n1=100;n2=100;s1=0.5;s2=2;ʱ�ĺ���ͼ��
x0=10;y0=10;
options =odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-5]);
[T,Y]=ode45('zhongqun',[0 50],[x0 y0],options);
grid on
axis equal
plot(T,Y(:,1),'b-',T,Y(:,2),'r-')
title('r1=1;r2=1;s1=0.5;s2=2;n1=100;n2=100;x0=10;y0=10;')
h = legend('x(t)','y(t)',2);
end


%�Զ�����Ⱥ����
function dy=zhongqun(t,y)
syms r1 r2 s1 s2 n1
%r��n���費ͬ�Ĳ���ʱ���в�ͬ�Ľ�
r1=1;r2=1;	
n1=100;n2=100;
s1=0.5;s2=2;
dy=zeros(2,1);
dy(1)=r1*y(1)*(1-y(1)/n1-s1*y(2)/n2);
dy(2)=r2*y(2)*(1-s2*y(1)/n1- y(2)/n2);
%ע�⣺�ڴ˺����У��ı�r1=;r2=;n1=;n2=;s1=;s2=;��ֵ���ﵽ���Ҫ��
end
