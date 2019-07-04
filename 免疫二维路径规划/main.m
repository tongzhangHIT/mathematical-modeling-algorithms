%% 清空环境
clc % 清屏
clear all; % 删除workplace变量
%close all; % 关掉显示图形窗口

%% 障碍物数据
position = load('barrier.txt');%改
plot([0,300],[0,300],'.');%改
hold on
B = load('barrier.txt');%改
xlabel('km','fontsize',12)
ylabel('km','fontsize',12)
title('二维规划空间','fontsize',12)
%% 描述起点和终点
S = [0,230];%改
T = [280,50];%改
plot([S(1),T(1)],[S(2),T(2)],'r.','MarkerSize',25);

% 图形标注
text(S(1)+2,S(2),'S');
text(T(1)+2,T(2),'T');
 
%% 描绘障碍物图形
fill(position(1:4,1),position(1:4,2),[0,1,0]);%改
fill(position(5:8,1),position(5:8,2),[0,1,0]);%改
%fill(position(9:12,1),position(9:12,2),[0,1,0]);%改
%fill(position(13:15,1),position(13:15,2),[0,1,0]);%改

% 下载链路端点数据
L = load('lines.txt');%改
 
%% 描绘线及中点
v = zeros(size(L));
for i=1:length(L)
    plot([position(L(i,1),1),position(L(i,2),1)],[position(L(i,1),2)...
        ,position(L(i,2),2)],'color','r','LineStyle','--','linewidth',2);
    v(i,:) = (position(L(i,1),:)+position(L(i,2),:))/2;
    plot(v(i,1),v(i,2),'*');
    text(v(i,1)+2,v(i,2),strcat('v',num2str(i)));
end
 
%% 描绘可行路径
sign = load('matrix.txt');%改
[n,m]=size(sign);
 
for i=1:n
    
    if i == 1
        for k=1:m-1
            if sign(i,k) == 1
                plot([S(1),v(k-1,1)],[S(2),v(k-1,2)],'color',...
                    'black','Linewidth',2,'LineStyle','-');
            end
        end
        continue;
    end
    
    for j=2:i
        if i == m
            if sign(i,j) == 1
                plot([T(1),v(j-1,1)],[T(2),v(j-1,2)],'color',...
                    'black','Linewidth',2,'LineStyle','-');
            end
        else
            if sign(i,j) == 1
                plot([v(i-1,1),v(j-1,1)],[v(i-1,2),v(j-1,2)],...
                    'color','black','Linewidth',2,'LineStyle','-');
            end
        end
    end
end
position1=ones(length(v)+2,2);
vv(1,:) = S;
vv(length(v),:) = T;
vv(2:length(v)+1,:) = v;

path = DijkstraPlan(position1,sign);
%path = DijkstraPlan(position,sign);
j = path(length(v)+2);
plot([T(1),v(j-1,1)],[T(2),v(j-1,2)],'color','g','LineWidth',3,'LineStyle','-.');
i = path(length(v)+2);
j = path(i);
count = 0;
while true
    plot([v(i-1,1),v(j-1,1)],[v(i-1,2),v(j-1,2)],'color','b','LineWidth',3,'LineStyle','-.');
    count = count + 1;
    i = j;
    j = path(i);
    if i == 1 || j==1
        break;
    end
end
plot([S(1),v(i-1,1)],[S(2),v(i-1,2)],'color','g','LineWidth',3,'LineStyle','-.');

path;
count = count+3;
pathtemp(count) = length(v)+2;
j =  length(v)+2;
for i=2:count
    pathtemp(count-i+1) = path(j);
    j = path(j);
end
path = pathtemp


%% 蚁群算法参数初始化
pathCount = length(path)-2;          %经过线段数量
pheCacuPara=2; %改                      %信息素计算参数
pheThres = 0.8; %改                     %信息素选择阈值
pheUpPara=[0.1 0.0003];  %改            %信息素更新参数
qfz= zeros(pathCount,10);            %启发值

phePara = ones(pathCount,10)*pheUpPara(2);         %信息素
qfzPara1 = ones(10,1)*0.5;           %启发信息参数
qfzPara2 = 1.1;                      %启发信息参数
m=10;    %改                            %种群数量
NC=500;  %改                            %循环次数
pathk = zeros(pathCount,m);          %搜索结果记录
shortestpath = zeros(1,NC);          %进化过程记录
 
%% 初始最短路径
dijpathlen = 0;
vv = zeros( length(v)+2,2);
vv(1,:) = S;
vv( length(v)+2,:) = T;
vv(2: length(v)+1,:) = v;
for i=1:pathCount-1
dijpathlen = dijpathlen + sqrt((vv(path(i),1)-vv(path(i+1),1))^2+(vv(path(i),2)-vv(path(i+1),2))^2);
end
LL = dijpathlen;
 
%% 经过的链接线
lines = zeros(pathCount,4);
for i = 1:pathCount
    lines(i,1:2) = B(L(path(i+1)-1,1),:);
    lines(i,3:4) = B(L(path(i+1)-1,2),:);
end
 
%% 循环搜索
for num = 1:NC
    
    %% 蚂蚁迭代寻优一次
    for i=1:pathCount
        for k=1:m
            q = rand();
            qfz(i,:) = (qfzPara2-abs((1:10)'/10-qfzPara1))/qfzPara2; %启发信息
            if q<=pheThres%选择信息素最大值
                arg = phePara(i,:).*(qfz(i,:).^pheCacuPara);
                j = find(arg == max(arg));
                pathk(i,k) = j(1);
            else  % 轮盘赌选择
                arg = phePara(i,:).*(qfz(i,:).^pheCacuPara);
                sumarg = sum(arg);
                qq = (q-pheThres)/(1-pheThres);
                qtemp = 0;
                j = 1;
                while qtemp < qq
                    qtemp = qtemp + (phePara(i,j)*(qfz(i,j)^pheCacuPara))/sumarg;
                    j=j+1;
                end
                j=j-1;
                pathk(i,k) = j(1);
            end
            % 信息素更新
            phePara(i,j) = (1-pheUpPara(1))*phePara(i,j)+pheUpPara(1)*pheUpPara(2);
        end
    end
    
    %% 计算路径长度
    len = zeros(1,k);
    for k=1:m
        Pstart = S;%起点
        Pend = lines(1,1:2) + (lines(1,3:4)-lines(1,1:2))*pathk(1,k)/10;
        for l=1:pathCount
            len(1,k) = len(1,k)+sqrt(sum((Pend-Pstart).^2));
            Pstart = Pend;
            if l<pathCount
                Pend = lines(l+1,1:2) + (lines(l+1,3:4)-lines(l+1,1:2))*pathk(l+1,k)/10;
            end
        end
        Pend = T;%终点
        len(1,k) = len(1,k)+sqrt(sum((Pend-Pstart).^2));
    end
    
    %% 更新信息素
    % 寻找最短路径
    minlen = min(len);
    minlen = minlen(1);
    minant = find(len == minlen);
    minant = minant(1);
    
    % 更新全局最短路径
    if minlen < LL
        LL = minlen;
    end
    
    % 更新信息素
    for i=1:pathCount
        phePara(i,pathk(i,minant)) = (1-pheUpPara(1))* phePara(i,pathk(i,minant))+pheUpPara(1)*(1/minlen);
    end
    shortestpath(num) = minlen;
end

figure;
plot(1:NC,shortestpath,'linewidth',2);
grid on
ylabel('路径总长度');
xlabel('迭代次数');
min(shortestpath)
