%----------------------------------------
%   基于贝叶斯判别的机场航班延误因素分析
%----------------------------------------

clc,clear,close all
load('sourcedata.mat');%加载用到的数据表格,原数据来源，与编程无关
load data.mat
load('datatest.mat');
n=size(data);

% 创建朴素贝叶斯分类器对象ObjBayes

training=data(1:103,1:5);%训练对象--可更改
group=data(1:103,6);%结果--可更改
ObjBayes = NaiveBayes.fit(training,group,'Distribution','kernel');%根据数据分类创建贝叶斯预测集


% 利用所创建的朴素贝叶斯分类器对象ObjBayes，对训练样本进行判别

pre0 = ObjBayes.predict(training);%根据training进行结果的预测-用处不大
disp '贝叶斯分类器训练数据和实际结果是否相等，相等为1，否则为0'%用处不大
isequal(pre0, group)  % 判断判别结果pre0与分组向量group是否相等—用处不大

pre1 = ObjBayes.predict(data(1:103,1:5));

% isequal(pre1, data(71:103,6))  % 判断判别结果pre0与分组向量group是否相等
figure,
subplot(211),bar(data(:,6));figure(gcf);axis tight,box off,grid on
title('原始数据---> 用于训练网络---103组数据 ---实际延误率')
subplot(212),bar(pre1);figure(gcf);axis tight,box off,grid on
title('贝叶斯网络训练结果---预测延误率')




%%
% 对于样本进行预测
test=datatest(:,1:5);
datatestresult=datatest(:,6);
pre2 = ObjBayes.predict(test);
figure,
%isequal(pre1, datatestresult)  % 判断判别结果pre0与分组向量group是否相等
subplot(211),bar(datatest(:,6));figure(gcf);axis tight,box off,grid on
title('输入待检验的数据，实际结果')
subplot(212),bar(pre2);figure(gcf);axis tight,box off,grid on
title('贝叶斯网络训练结果')



