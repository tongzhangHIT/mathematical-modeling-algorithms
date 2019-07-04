clear;
close all;
clc;

x=-0.5:0.01:0.5;  %-----------可替换
y=-0.5:0.01:0.5;  %-----------可替换
for i=1:101
    for j=1:101
z(i,j)=-20*exp(-0.2*sqrt((x(i)^2+y(j)^2)/2))-exp((cos(2*pi*x(i))+cos(2*pi*y(j)))/2)+20+2.71928;  %-----------可替换
    end
end
mesh(z);
hold on;
n=100;   %粒子群粒子个数  %-----------可替换

%初始化粒子群，定义结构体
%结构体中八个元素，分别是粒子坐标，粒子速度，粒子适应度，粒子最佳适应度，粒子最佳坐标
par=struct([]);
for i=1:n
    par(i).x=-0.5+rand();   %[-100 100]对x位置随机初始化  %-----------可替换
    par(i).y=-0.5+rand();   %[-100 100]对y位置随机初始化  %-----------可替换
    par(i).vx=-1+2*rand();      %[-1 1]对vx速度随机初始化
    par(i).vy=-1+2*rand();      %[-1 1]对vy速度随机初始化
    par(i).fit=0;               %粒子适应度为0初始化
    par(i).bestfit=0;           %粒子最佳适应度为0初始化
    par(i).bestx=par(i).x;      %粒子x最佳位置初始化
    par(i).besty=par(i).y;      %粒子y最佳位置初始化
end
par_best=par(1);    %初始化粒子群中最佳粒子

for k=1:20    
    plot3(par_best.x,par_best.y,par_best.fit,'g*'); %画出最佳粒子的位置  %-----------可替换
    for p=1:n
        [par(p),par_best]=update_par(par(p),par_best);  %更新每个粒子信息         
    end  
end