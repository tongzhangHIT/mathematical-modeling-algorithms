clear all;
load genezhongxin.mat;load gene.mat;

P=genezhongxin;
P(3,:)=gene(50,:);
[Pn,Ps]=mapminmax(P);
P=Pn;
P=P';

data=gene;
T=data(41:60,:);
[Tn,Ts]=mapminmax(T);
%% 网络建立和训练
% newsom建立SOM网络。minmax（P）取输入的最大最小值。竞争层为6*6=36个神经元
net=newsom(minmax(P),[5 5]);%神经元个数
plotsom(net.layers{1}.positions)
% 5次训练的步数
a=[1000];%训练次数
% 随机初始化一个1*10向量。
yc=rands(1,3);%7种训练次数，8种结果
%% 进行训练
% 训练次数为10次
net.trainparam.epochs=a(1);
% 训练网络和查看分类
net=train(net,P);
y=sim(net,P);
yc(1,:)=vec2ind(y);
%figure
plotsom(net.IW{1,1},net.layers{1}.distances)

t=T';
%t=T(1,:)';
r=sim(net,t);
% 变换函数 将单值向量转变成下标向量。
rr=vec2ind(r)


figure,plotsomtop(net)
% 查看临近神经元直接的距离情况
figure,plotsomnd(net)
% 查看每个神经元的分类情况
figure,plotsomhits(net,P)
