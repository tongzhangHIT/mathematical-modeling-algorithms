clc;
clear all;
close all;

f=@(x)x.*sin(x)+x.*cos(2.*x);%------可替换

lb=0;%------可替换
ub=10;%------可替换

%设定粒子的初始位置
x0=[0 1 3 6 8 10];%------可替换

hf=figure;
for i=1:6
    %fmincon 有约束的非线性最小化
    %fmincon(fun,x,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
    x(i)=fmincon(f,x0(i),[],[],[],[],lb,ub,[],...
        optimset('Algorithm','SQP','Disp','none'));
    subplot(2,3,i)
    ezplot(f,[lb ub]);
    hold on;
    plot(x0(i),f(x0(i)),'k+');
    plot(x(i),f(x(i)),'ro');
    hold off
    title(['Starting at',num2str(x0(i))]);%绘制起始点图像
    if i==1||i==4
        ylabel('xsin(x)+xcos(x)');
    end
end

%方案2--fminbnd 单变量边界非线性函数最小值求解
x2=fminbnd(f,lb,ub)
figure
ezplot(f,[lb ub]);
hold on
plot(x2,f(x2),'ro');
hold off
ylabel('xsin(x)+xcos(x)');
title('Solution using fminbnd');

