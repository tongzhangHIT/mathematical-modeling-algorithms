%% SOM神经网络的数据分类--柴油机故障诊断
clc % 清屏
clear all; % 删除workplace变量
close all; % 关掉显示图形窗口
format short
% Initial
%% 录入输入数据
% 载入数据
load p;
P=P(:,1:3);

%转置后符合神经网络的输入格式
P=P';

%% 网络建立和训练
% newsom建立SOM网络。minmax（P）取输入的最大最小值。竞争层为6*6=36个神经元
net=newsom(minmax(P),[5 5]);%神经元个数
plotsom(net.layers{1}.positions)
% 5次训练的步数
a=[1000 ];%训练次数
% 随机初始化一个1*10向量。
yc=rands(7,3);%7种训练次数，8种结果
%% 进行训练
% 训练次数为10次
net.trainparam.epochs=a(1);
% 训练网络和查看分类
net=train(net,P);
y=sim(net,P);
yc(1,:)=vec2ind(y);
%figure
plotsom(net.IW{1,1},net.layers{1}.distances)





%% 网络神经元分布情况
% 查看网络拓扑学结构
figure,plotsomtop(net)
% 查看临近神经元直接的距离情况
figure,plotsomnd(net)
% 查看每个神经元的分类情况
figure,plotsomhits(net,P)