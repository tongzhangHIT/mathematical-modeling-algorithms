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
%% ���罨����ѵ��
% newsom����SOM���硣minmax��P��ȡ����������Сֵ��������Ϊ6*6=36����Ԫ
net=newsom(minmax(P),[5 5]);%��Ԫ����
plotsom(net.layers{1}.positions)
% 5��ѵ���Ĳ���
a=[1000];%ѵ������
% �����ʼ��һ��1*10������
yc=rands(1,3);%7��ѵ��������8�ֽ��
%% ����ѵ��
% ѵ������Ϊ10��
net.trainparam.epochs=a(1);
% ѵ������Ͳ鿴����
net=train(net,P);
y=sim(net,P);
yc(1,:)=vec2ind(y);
%figure
plotsom(net.IW{1,1},net.layers{1}.distances)

t=T';
%t=T(1,:)';
r=sim(net,t);
% �任���� ����ֵ����ת����±�������
rr=vec2ind(r)


figure,plotsomtop(net)
% �鿴�ٽ���Ԫֱ�ӵľ������
figure,plotsomnd(net)
% �鿴ÿ����Ԫ�ķ������
figure,plotsomhits(net,P)
