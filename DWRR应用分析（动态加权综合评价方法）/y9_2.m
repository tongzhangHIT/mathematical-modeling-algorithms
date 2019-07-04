clc% 清屏
clear all;%删除workplace变量
close all;%关掉显示图形窗口
load('A.mat')
ca=[0.1,0.335,0.1];%改%中值
ca1=[0.5341  0.4492  0.5341];%改%δ的值
B=[0.2 0.6 1;0.67 0.67 1;0.2 0.6 1];%改%标准分级
a=size(A,1);%A的行数
b=size(A,2);%A的列数
MX=max(A);%A的每列的最大值
MN=min(A);%A的每列的最小值
f1=1;%B的行变量
f2=1;%B的列变量
for j=1:b%将A里面的值标准化
    for i=1:a
        A(i,j)=(A(i,j)-MN(j))/(MX(j)-MN(j));
    end
end
X=zeros(a,4);%改（4个城市）%建立放置每个城市每天空气质量的值矩阵
for i=1:a%计算每个城市每天的空气质量
    h=1;
    k=1;
    flag0=1;%算完每个城市每天的空气质量后，跳到下一个城市（h=h+1）
    for j=1:b 
        for flag=1:3%改（3种参数）
            if A(i,j)<=ca(k)
              w=0;
              X(i,h)=X(i,h);
              break;
            elseif A(i,j)>ca(k)&&A(i,j)<B(f1,f2)
                w=1-exp(-(A(i,j)-ca(k)/ca1(k))^2);
                X(i,h)=X(i,h)+w*A(i,j);%加权
                break;
            else
                k=k+1;%对应每个中值
                f2=f2+1;%对应每个等级值
            end
        end
        f1=f1+1;%对应每个属性
        k=1;%每循环完一次，跳到第一个中值
        f2=1;%每循环完一次，又从第一级开始
        if flag0==3%改（3种参数）%判断是否算完了一个城市每天的空气质量
            h=h+1;
            flag0=1;
            f1=1;
        else
        flag0=flag0+1;
        end
    end
end
a=size(X,1);%X的行数
b=size(X,2);%X的列数
h=1;
C=zeros(1,b);%建立总体排序矩阵
for k=1:b%将所有城市总体排序
    for i=1:a
        for j=1:b
            if(A(h)>A(i,j))
                C(k)= C(k)+1;
            else
            end
        end
        h=h+1;
    end
end
C % Borda数
X;