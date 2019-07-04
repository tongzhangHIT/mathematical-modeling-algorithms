clc;
clear ;
close all;
global E  F  B  K 
E=0:0.1:2;%----�ɽ��и��� r1�ı仯��Χ
E=E';
F=0:0.1:2;%----�ɽ��и��� r2�ı仯��Χ
F=F';
B=1;K=1;a=1;b=1;c=1;d=1;
S=zeros(140,2);         %��¼�������볤�뾨������ȶ�״̬����Ϊ���s1,s2ֵ;
H=zeros(250,2);         %��¼�����������������ֵ̬��Ϊx=150000,y=0��s1,s2ֵ;
U=zeros(250,2);         %���뾨�������������ֵ̬��Ϊx=0,y=400000��s1,s2ֵ;
Num=zeros(441,4);        %��¼�����볤�뾨������ȶ�״ֵ̬;
while B<22
     K=1;
     while K<22
        options =odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-5]);%----����ѡ�񣬿ɲ��Ķ�
        [T,Y] = ode45('zhongqun4',[0 2000],[5000 70000],options);%----�ɽ��и��� 
        [m,n]=size(Y);
        Num(a,1)=Y(m,1);
        if Num(a,1)<1
           Num(a,1)=0;
        end
        Num(a,2)=Y(m,2);
        if Num(a,2)<1
           Num(a,2)=0;
        end
       Num(a,3)=E(B);
       Num(a,4)=F(K);
        if ((Y(m,1)-1>0)&(Y(m,2)-1>0))==1
               S(b,1)=E(B);
               S(b,2)=F(K);
               b=b+1;
        end
         if (Y(m,2)-1)<0
               H(c,1)=E(B);
               H(c,2)=F(K);
               c=c+1;
         end
         if (Y(m,1)-1)<0
               U(d,1)=E(B);
               U(d,2)=F(K);
               d=d+1;
        end
        a=a+1;
        K=K+1;
     end
     B=B+1;
end
%--------------------------��ͼ-------------------------------------
[s1,s2]=meshgrid(Num(:,3),Num(:,4));
LANJING=griddata(Num(:,3),Num(:,4),Num(:,1),s1,s2,'v4');
CHANGXUJING=griddata(Num(:,3),Num(:,4),Num(:,2),s1,s2,'v4');
figure(1)
mesh(s1,s2,LANJING);
title('s1��s2�ı�ʱ�����������仯ͼ��');
figure(2)
mesh(s1,s2,CHANGXUJING);
title('s1��s2�ı�ʱ�����뾨�����仯ͼ��');


