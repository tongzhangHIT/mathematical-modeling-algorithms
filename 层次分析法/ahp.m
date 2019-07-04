clear all;clc;
a=[1 1 1 4 1 1/2
    1 1 2 4 1 1/2
    1 1/2 1 5 3 1/2
    1/4 1/4 1/5 1 1/3 1/3
    1 1 1/3 3 1 1
    2 2 2 3 3 1];%准则层相互因子权重矩阵
[x,y]=eig(a);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci1=(lamda-6)/5;%一致性指标：ci=（lambda-n）/(n-1)
cr1=ci1/1.24%一致性比例：cr=ci/ri.要求cr<0.1
 %不同n值对应不同ri
  %1  2  3    4   5    6    7    8     9
  %0  0 0.58 0.9 1.12 1.24 1.32 1.41 1.45
w1=x(:,1)/sum(x(:,1))%准则层权重结果
%准则层1
b1=[1 1/4 1/2
    4 1 3
    2 1/3 1];%方案层对准则层1相互因子权重矩阵
[x,y]=eig(b1);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci21=(lamda-3)/2;%ci=（lambda-n）/(n-1)
ri2=0.58;
cr21=ci21/ri2%cr=ci/ri
w21=x(:,1)/sum(x(:,1))%准则层1权值
%准则层2
b2=[1 1/4 1/5
    4 1 1/2
    5 2 1];%方案层对准则层2相互因子权重矩阵
[x,y]=eig(b2);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci22=(lamda-3)/2;%ci=（lambda-n）/(n-1)
cr22=ci22/ri2%cr=ci/ri
w22=x(:,1)/sum(x(:,1))%准则层2权值
%准则层3
b3=[1 3 1/3
    1/3 1 1/7
    3 7 1];%方案层对准则层3相互因子权重矩阵
[x,y]=eig(b3);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci23=(lamda-3)/2;%ci=（lambda-n）/(n-1)
cr23=ci23/ri2%cr=ci/ri
w23=x(:,1)/sum(x(:,1))%准则层3权值
%准则层4
b4=[1 1/3 5
    3 1 7
    1/5 1/7 1];%方案层对准则层4相互因子权重矩阵
[x,y]=eig(b4);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci24=(lamda-3)/2;%ci=（lambda-n）/(n-1)
cr24=ci24/ri2%cr=ci/ri
w24=x(:,1)/sum(x(:,1))%准则层4权值
%准则层5
b5=[1 1 7
    1 1 7
    1/7 1/7 1];%方案层对准则层5相互因子权重矩阵
[x,y]=eig(b5);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci25=(lamda-3)/2;%ci=（lambda-n）/(n-1)
cr25=ci25/ri2%cr=ci/ri
w25=x(:,1)/sum(x(:,1))%准则层5权值
%准则层6
b6=[1 7 9
    1/7 1 1
    1/9 1 1];%方案层对准则层6相互因子权重矩阵
[x,y]=eig(b6);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci26=(lamda-3)/2;%ci=（lambda-n）/(n-1)
cr26=ci21/ri2%cr=ci/ri
w26=x(:,1)/sum(x(:,1))%准则层6权值
%总排序
w_sum=[w21,w22,w23,w24,w25,w26]*w1%方案层总排序权值
ci=[ci21,ci22,ci23,ci24,ci25,ci26];
cr=ci*w1/sum(ri2*w1)%总排序一致性比例

