clc% ����
clear all;%ɾ��workplace����
close all;%�ص���ʾͼ�δ���
load('A.mat')
ca=[0.1,0.335,0.1];%��%��ֵ
ca1=[0.5341  0.4492  0.5341];%��%�ĵ�ֵ
B=[0.2 0.6 1;0.67 0.67 1;0.2 0.6 1];%��%��׼�ּ�
a=size(A,1);%A������
b=size(A,2);%A������
MX=max(A);%A��ÿ�е����ֵ
MN=min(A);%A��ÿ�е���Сֵ
f1=1;%B���б���
f2=1;%B���б���
for j=1:b%��A�����ֵ��׼��
    for i=1:a
        A(i,j)=(A(i,j)-MN(j))/(MX(j)-MN(j));
    end
end
X=zeros(a,4);%�ģ�4�����У�%��������ÿ������ÿ�����������ֵ����
for i=1:a%����ÿ������ÿ��Ŀ�������
    h=1;
    k=1;
    flag0=1;%����ÿ������ÿ��Ŀ���������������һ�����У�h=h+1��
    for j=1:b 
        for flag=1:3%�ģ�3�ֲ�����
            if A(i,j)<=ca(k)
              w=0;
              X(i,h)=X(i,h);
              break;
            elseif A(i,j)>ca(k)&&A(i,j)<B(f1,f2)
                w=1-exp(-(A(i,j)-ca(k)/ca1(k))^2);
                X(i,h)=X(i,h)+w*A(i,j);%��Ȩ
                break;
            else
                k=k+1;%��Ӧÿ����ֵ
                f2=f2+1;%��Ӧÿ���ȼ�ֵ
            end
        end
        f1=f1+1;%��Ӧÿ������
        k=1;%ÿѭ����һ�Σ�������һ����ֵ
        f2=1;%ÿѭ����һ�Σ��ִӵ�һ����ʼ
        if flag0==3%�ģ�3�ֲ�����%�ж��Ƿ�������һ������ÿ��Ŀ�������
            h=h+1;
            flag0=1;
            f1=1;
        else
        flag0=flag0+1;
        end
    end
end
a=size(X,1);%X������
b=size(X,2);%X������
h=1;
C=zeros(1,b);%���������������
for k=1:b%�����г�����������
    for i=1:a
        for j=1:b
            if(A(h)>A(i,j))
                C(k)= C(k)+1;
            else
            end
        end
        h=h+1;
    end
end
C % Borda��
X;