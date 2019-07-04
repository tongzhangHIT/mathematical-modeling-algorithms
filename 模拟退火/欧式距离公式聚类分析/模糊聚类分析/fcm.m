clear all
close all
clc

%Initial
data=rand(100,2);
[center,U,obj_fcn]=fcm(data,2);
plot(data(:,1),data(:,2),'o');
maxU=max(U);
index1=find(U(1,:)==maxU);
index2=find(U(2,:)==maxU);
line(data(index1,1),data(index1,2),'linestyle','none','marker','o','color','g');
line(data(index2,1),data(index2,2),'linestyle','none','marker','o','color','r');