clc;
clear;
close all
%各类食品的价格预测
% x由workplace载入

load('x.mat')
x=x';

for i=1:37                      %42类38组数据进行一次移动
      for j=1:42
       m1(i+2,j)=(x(i,j)+x(i+1,j)+x(i+2,j))/3;
end
end

for i=3:37                      %进行二次移动
     for j=1:42
m2(i+2,j)=(m1(i,j)+m1(i+1,j)+m1(i+2,j))/3;
     end
end

a=2*m1-m2;                      %参数有公式确定
b=m1-m2;

for i=1:42
      for j=1:6
          y(i,j)=a(39,i)+b(39,i)*j;      %完成数据的预测
      end
end

