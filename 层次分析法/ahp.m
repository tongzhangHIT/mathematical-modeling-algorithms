clear all;clc;
a=[1 1 1 4 1 1/2
    1 1 2 4 1 1/2
    1 1/2 1 5 3 1/2
    1/4 1/4 1/5 1 1/3 1/3
    1 1 1/3 3 1 1
    2 2 2 3 3 1];%׼����໥����Ȩ�ؾ���
[x,y]=eig(a);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci1=(lamda-6)/5;%һ����ָ�꣺ci=��lambda-n��/(n-1)
cr1=ci1/1.24%һ���Ա�����cr=ci/ri.Ҫ��cr<0.1
 %��ͬnֵ��Ӧ��ͬri
  %1  2  3    4   5    6    7    8     9
  %0  0 0.58 0.9 1.12 1.24 1.32 1.41 1.45
w1=x(:,1)/sum(x(:,1))%׼���Ȩ�ؽ��
%׼���1
b1=[1 1/4 1/2
    4 1 3
    2 1/3 1];%�������׼���1�໥����Ȩ�ؾ���
[x,y]=eig(b1);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci21=(lamda-3)/2;%ci=��lambda-n��/(n-1)
ri2=0.58;
cr21=ci21/ri2%cr=ci/ri
w21=x(:,1)/sum(x(:,1))%׼���1Ȩֵ
%׼���2
b2=[1 1/4 1/5
    4 1 1/2
    5 2 1];%�������׼���2�໥����Ȩ�ؾ���
[x,y]=eig(b2);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci22=(lamda-3)/2;%ci=��lambda-n��/(n-1)
cr22=ci22/ri2%cr=ci/ri
w22=x(:,1)/sum(x(:,1))%׼���2Ȩֵ
%׼���3
b3=[1 3 1/3
    1/3 1 1/7
    3 7 1];%�������׼���3�໥����Ȩ�ؾ���
[x,y]=eig(b3);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci23=(lamda-3)/2;%ci=��lambda-n��/(n-1)
cr23=ci23/ri2%cr=ci/ri
w23=x(:,1)/sum(x(:,1))%׼���3Ȩֵ
%׼���4
b4=[1 1/3 5
    3 1 7
    1/5 1/7 1];%�������׼���4�໥����Ȩ�ؾ���
[x,y]=eig(b4);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci24=(lamda-3)/2;%ci=��lambda-n��/(n-1)
cr24=ci24/ri2%cr=ci/ri
w24=x(:,1)/sum(x(:,1))%׼���4Ȩֵ
%׼���5
b5=[1 1 7
    1 1 7
    1/7 1/7 1];%�������׼���5�໥����Ȩ�ؾ���
[x,y]=eig(b5);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci25=(lamda-3)/2;%ci=��lambda-n��/(n-1)
cr25=ci25/ri2%cr=ci/ri
w25=x(:,1)/sum(x(:,1))%׼���5Ȩֵ
%׼���6
b6=[1 7 9
    1/7 1 1
    1/9 1 1];%�������׼���6�໥����Ȩ�ؾ���
[x,y]=eig(b6);
eigenvalue=diag(y);
lamda=eigenvalue(1);
ci26=(lamda-3)/2;%ci=��lambda-n��/(n-1)
cr26=ci21/ri2%cr=ci/ri
w26=x(:,1)/sum(x(:,1))%׼���6Ȩֵ
%������
w_sum=[w21,w22,w23,w24,w25,w26]*w1%������������Ȩֵ
ci=[ci21,ci22,ci23,ci24,ci25,ci26];
cr=ci*w1/sum(ri2*w1)%������һ���Ա���

