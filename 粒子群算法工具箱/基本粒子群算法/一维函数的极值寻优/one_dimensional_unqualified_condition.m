clc;
clear all;
close all;

f=@(x)x.*sin(x)+x.*cos(2.*x);%------���滻

lb=0;%------���滻
ub=10;%------���滻

%�趨���ӵĳ�ʼλ��
x0=[0 1 3 6 8 10];%------���滻

hf=figure;
for i=1:6
    %fmincon ��Լ���ķ�������С��
    %fmincon(fun,x,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
    x(i)=fmincon(f,x0(i),[],[],[],[],lb,ub,[],...
        optimset('Algorithm','SQP','Disp','none'));
    subplot(2,3,i)
    ezplot(f,[lb ub]);
    hold on;
    plot(x0(i),f(x0(i)),'k+');
    plot(x(i),f(x(i)),'ro');
    hold off
    title(['Starting at',num2str(x0(i))]);%������ʼ��ͼ��
    if i==1||i==4
        ylabel('xsin(x)+xcos(x)');
    end
end

%����2--fminbnd �������߽�����Ժ�����Сֵ���
x2=fminbnd(f,lb,ub)
figure
ezplot(f,[lb ub]);
hold on
plot(x2,f(x2),'ro');
hold off
ylabel('xsin(x)+xcos(x)');
title('Solution using fminbnd');

