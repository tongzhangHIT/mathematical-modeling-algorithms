clc;
clear;
close all
load('jlx.mat');
x=jlx;
BX=zscore(x);%��׼�����ݾ���
Y=pdist(x);%��ŷ�Ͼ����������֮��ľ��룻
D=squareform(Y);%ŷ�Ͼ������
Z=linkage(Y);%��̾��뷨��
% T=cluster(Z,4);
[H,T]=dendrogram(Z,'colorthreshold','default');