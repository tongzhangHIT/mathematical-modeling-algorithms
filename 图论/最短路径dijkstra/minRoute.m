function [S,D]=minRoute(i,m,W,opt)
% ͼ���������������·����Dijkstra�㷨M����
% ��ʽ[S,D]=minRoute(i,m,W,opt)
% iΪ���·������ʼ�㣬mΪͼ��������WΪͼ�Ĵ�Ȩ�ڽӾ��󣬲�����
% �ߵ�������֮���Ȩ��inf��ʾ��S��ÿһ�д��ϵ��¼�¼�˴�ʼ�㵽�յ�
% �����·�������������š�opt=0(ȱʡֵ)ʱ��S���յ���Ŵ�С������
% ʾ�����opt=1ʱ��S�����·����С������ʾ�����D��һ����������Ӧ
% ��¼��S������ʾ·���Ĵ�С��

%����ʾ��
%n=6;
%W=inf*ones(6);
%W(1,[2,3,4,5,6])=[25,32,44,62,87];
%W(2,[3,4,5,6])=[25,32,44,62];
%W(3,[4,5,6])=[27,34,46];
%W(4,[5,6])=[27,34];
%W(5,6)=28;
if nargin<4
    opt=0;
end
dd=[];tt=[];
ss=[];ss(1,1)=i;
V=1:m;V(i)=[];
dd=[0;i];
kk=2;
[mdd,ndd]=size(dd);
while ~isempty(V)
    [tmpd,j]=min(W(i,V));
    tmpj=V(j);
    for k=2:ndd
        [tmp1,jj]=min(dd(1,k)+W(dd(2,k),V));
        tmp2=V(jj);
        tt(k-1,:)=[tmp1,tmp2,jj];
    end
    tmp=[tmpd,tmpj,j;tt];
    [tmp3,tmp4]=min(tmp(:,1));
    if tmp3==tmpd
        ss(1:2,kk)=[i;tmp(tmp4,2)];
    else
        tmp5=find(ss(:,tmp4)~=0);
        tmp6=length(tmp5);
        if dd(2,tmp4)==ss(tmp6,tmp4)
            ss(1:tmp6+1,kk)=[ss(tmp5,tmp4);tmp(tmp4,2)];
        else
            ss(1:3,kk)=[i;dd(2,tmp4);tmp(tmp4,2)];
        end
    end
    dd=[dd,[tmp3;tmp(tmp4,2)]];
    V(tmp(tmp4,3))=[];
    [mdd,ndd]=size(dd);
    kk=kk+1;
end
if opt==1
    [tmp,t]=sort(dd(2,:));
    S=ss(:,t);
    D=dd(1,t);
else
    S=ss;
    D=dd(1,:);
end
