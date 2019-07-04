%----------------------------------------
%   ���ڱ�Ҷ˹�б�Ļ��������������ط���
%----------------------------------------

clc,clear,close all
load('sourcedata.mat');%�����õ������ݱ��,ԭ������Դ�������޹�
load data.mat
load('datatest.mat');
n=size(data);

% �������ر�Ҷ˹����������ObjBayes

training=data(1:103,1:5);%ѵ������--�ɸ���
group=data(1:103,6);%���--�ɸ���
ObjBayes = NaiveBayes.fit(training,group,'Distribution','kernel');%�������ݷ��ഴ����Ҷ˹Ԥ�⼯


% ���������������ر�Ҷ˹����������ObjBayes����ѵ�����������б�

pre0 = ObjBayes.predict(training);%����training���н����Ԥ��-�ô�����
disp '��Ҷ˹������ѵ�����ݺ�ʵ�ʽ���Ƿ���ȣ����Ϊ1������Ϊ0'%�ô�����
isequal(pre0, group)  % �ж��б���pre0���������group�Ƿ���ȡ��ô�����

pre1 = ObjBayes.predict(data(1:103,1:5));

% isequal(pre1, data(71:103,6))  % �ж��б���pre0���������group�Ƿ����
figure,
subplot(211),bar(data(:,6));figure(gcf);axis tight,box off,grid on
title('ԭʼ����---> ����ѵ������---103������ ---ʵ��������')
subplot(212),bar(pre1);figure(gcf);axis tight,box off,grid on
title('��Ҷ˹����ѵ�����---Ԥ��������')




%%
% ������������Ԥ��
test=datatest(:,1:5);
datatestresult=datatest(:,6);
pre2 = ObjBayes.predict(test);
figure,
%isequal(pre1, datatestresult)  % �ж��б���pre0���������group�Ƿ����
subplot(211),bar(datatest(:,6));figure(gcf);axis tight,box off,grid on
title('�������������ݣ�ʵ�ʽ��')
subplot(212),bar(pre2);figure(gcf);axis tight,box off,grid on
title('��Ҷ˹����ѵ�����')



