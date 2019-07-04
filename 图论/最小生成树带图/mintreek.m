function [wt,pp]=mintreek(n,W)
% 图论中最小生成树 Kruskal算法及画图程序M函数
% 格式[wt,pp]=mintreek(n,W):n为图顶点数,W为图的带权邻接矩阵,
% 不构成边的两顶点之间的权用Inf表示.显示最小生成树的边及顶点, 
% wt为最小生成树的权,pp(:,1,2)为最小生成树边的两顶点,pp(:,3)
% 为最小生成树的边权,pp(:,4)为最小生成树边的序号;

%输入示例
%n=5;
%W=inf*ones(5);
%W(1,[2,3,4])=[1,7,3];
%W(2,[3,5])=[6,4];
%W(3,[4,5])=[8,5];
%W(4,5)=2;
tmpa=find(W~=inf);
[tmpb,tmpc]=find(W~=inf);
w=W(tmpa);   % w是W中非inf元素按列构成的向量
e=[tmpb,tmpc];  % e的每一行元素表示一条边的两个顶点的序号
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
for i=1:length(p(:,3));   %显示顶点vi与边ej
    disp([' ','e',num2str(p(i,4)),'','(v',num2str(p(i,1)),'','v',num2str(p(i,2)),')']);
end
% 以下是画图程序
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
for i=1:nE  %画两两之间连线
    plot(x(1,e(i,:)),y(1,e(i,:)),'b');
end;
for i=1:length(p(:,4))%画最小树
    plot(x(1,pp(i,1:2)),y(1,pp(i,1:2)),'r');
end
text(-0.35,-1.2,['最小生成树的权为','',num2str(wt)]);
title('红色连线为最小生成树');
axis('off');
hold off;
