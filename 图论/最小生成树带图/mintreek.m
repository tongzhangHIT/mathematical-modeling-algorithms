function [wt,pp]=mintreek(n,W)
% ͼ������С������ Kruskal�㷨����ͼ����M����
% ��ʽ[wt,pp]=mintreek(n,W):nΪͼ������,WΪͼ�Ĵ�Ȩ�ڽӾ���,
% �����ɱߵ�������֮���Ȩ��Inf��ʾ.��ʾ��С�������ı߼�����, 
% wtΪ��С��������Ȩ,pp(:,1,2)Ϊ��С�������ߵ�������,pp(:,3)
% Ϊ��С�������ı�Ȩ,pp(:,4)Ϊ��С�������ߵ����;

%����ʾ��
%n=5;
%W=inf*ones(5);
%W(1,[2,3,4])=[1,7,3];
%W(2,[3,5])=[6,4];
%W(3,[4,5])=[8,5];
%W(4,5)=2;
tmpa=find(W~=inf);
[tmpb,tmpc]=find(W~=inf);
w=W(tmpa);   % w��W�з�infԪ�ذ��й��ɵ�����
e=[tmpb,tmpc];  % e��ÿһ��Ԫ�ر�ʾһ���ߵ�������������
[wa,wb]=sort(w);
E=[e(wb,:),wa,wb];
[nE,mE]=size(E);
temp=find(E(:,1)-E(:,2));
E=E(temp,:);
p=E(1,:);
k=length(E(:,1));
while (rank(E)>0)
    temp1=max(E(1,2),E(1,1));
    temp2=min(E(1,2),E(1,1));
    for i=1:k;
        if (E(i,1)==temp1)
            E(i,1)=temp2;
        end
        if (E(i,2)==temp1)
            E(i,2)=temp2;
        end
    end
    a=find(E(:,1)-E(:,2));
    E=E(a,:);
    if (rank(E)>0)
        p=[p;E(1,:)];
        k=length(E(:,1));
    end
end
wt=sum(p(:,3));
pp=[e(p(:,4),:),p(:,3:4)];
for i=1:length(p(:,3));   %��ʾ����vi���ej
    disp([' ','e',num2str(p(i,4)),'','(v',num2str(p(i,1)),'','v',num2str(p(i,2)),')']);
end
% �����ǻ�ͼ����
axis equal;
hold on
[x,y]=cylinder(1,n);
xm=min(x(1,:));
ym=min(y(1,:));
xx=max(x(1,:));
yy=max(y(1,:));
%axis([xm -abs(xm)*0.15,xx+abs(xx)*0.15,ym-abs(ym)*0.15,yy+abs(yy)*0.15]);
plot(x(1,:),y(1,:),'ko');
for i=1:n
    temp=[' v',int2str(i)];
    text(x(1,i),y(1,i),temp);
end
for i=1:nE  %������֮������
    plot(x(1,e(i,:)),y(1,e(i,:)),'b');
end;
for i=1:length(p(:,4))%����С��
    plot(x(1,pp(i,1:2)),y(1,pp(i,1:2)),'r');
end
text(-0.35,-1.2,['��С��������ȨΪ','',num2str(wt)]);
title('��ɫ����Ϊ��С������');
axis('off');
hold off;
