clc,clear
load('yx.mat') %ԭʼ�������������ķ�ʽ����ڴ��ı��ļ���
yt=yx; n=length(yt);
alpha=0.9; st1(1)=yt(1); st2(1)=yt(1);
for i=2:n
    st1(i)=alpha*yt(i)+(1-alpha)*st1(i-1);%����ָ��Ԥ��
    st2(i)=alpha*st1(i)+(1-alpha)*st2(i-1);
end
a=2*st1-st2;
b=alpha/(1-alpha)*(st1-st2);       %�����ɹ�ʽ����ó����޸���
yhat=a+b;
yhat=yhat';
str=char(['C',int2str(n+2)]);

%% ����ָ��ƽ��
clc,clear
load('yx.mat') %ԭʼ�������������ķ�ʽ����ڴ��ı��ļ���
yt=yx; n=length(yt);
alpha=0.9; st1_0=mean(yt(1:3)); st2_0=st1_0; st3_0=st1_0;
st1(1)=alpha*yt(1)+(1-alpha)*st1_0;
st2(1)=alpha*st1(1)+(1-alpha)*st2_0;
st3(1)=alpha*st2(1)+(1-alpha)*st3_0;
for i=2:n
    st1(i)=alpha*yt(i)+(1-alpha)*st1(i-1);  %��������ָ��Ԥ��
    st2(i)=alpha*st1(i)+(1-alpha)*st2(i-1);
    st3(i)=alpha*st2(i)+(1-alpha)*st3(i-1);
end
st1=[st1_0,st1]; 
st2=[st2_0,st2];
st3=[st3_0,st3];
a=3*st1-3*st2+st3;
b=0.5*alpha/(1-alpha)^2*((6-5*alpha)*st1-2*(5-4*alpha)*st2+(4-3*alpha)*st3);
c=0.5*alpha^2/(1-alpha)^2*(st1-2*st2+st3);  %����ֵ�ɹ�ʽȷ��
yhat=a+b+c;
yhat=yhat';
