%% SOM����������ݷ���--���ͻ��������
clc % ����
clear all; % ɾ��workplace����
close all; % �ص���ʾͼ�δ���
format short
% Initial
%% ¼����������
% ��������
load p;
P=P(:,1:3);

%ת�ú����������������ʽ
P=P';

%% ���罨����ѵ��
% newsom����SOM���硣minmax��P��ȡ����������Сֵ��������Ϊ6*6=36����Ԫ
net=newsom(minmax(P),[5 5]);%��Ԫ����
plotsom(net.layers{1}.positions)
% 5��ѵ���Ĳ���
a=[1000 ];%ѵ������
% �����ʼ��һ��1*10������
yc=rands(7,3);%7��ѵ��������8�ֽ��
%% ����ѵ��
% ѵ������Ϊ10��
net.trainparam.epochs=a(1);
% ѵ������Ͳ鿴����
net=train(net,P);
y=sim(net,P);
yc(1,:)=vec2ind(y);
%figure
plotsom(net.IW{1,1},net.layers{1}.distances)





%% ������Ԫ�ֲ����
% �鿴��������ѧ�ṹ
figure,plotsomtop(net)
% �鿴�ٽ���Ԫֱ�ӵľ������
figure,plotsomnd(net)
% �鿴ÿ����Ԫ�ķ������
figure,plotsomhits(net,P)