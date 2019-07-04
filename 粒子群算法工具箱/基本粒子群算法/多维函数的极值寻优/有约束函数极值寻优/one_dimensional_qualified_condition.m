clc
clear 
close all
warning off

%%参数的初始化
c1=1.49445;%---------可进行更改-------
c2=1.49445;%---------可进行更改-------
maxgen=200; %进化次数%---------可进行更改-------
sizepop=200;%种群规模%---------可进行更改-------

Vmax=1;     %粒子更新速度
Vmin=-1;

popmax=50;%种群
popmin=-50;

par_num=7;%维度

%%产生初始粒子和速度
for i=1:sizepop
    pop(i,:)=1.*rands(1,par_num);%初始化种群
    V(i,:)=1.*rands(1,par_num);%初始化速度
    %计算适应度
    fitness(i)=fun(pop(i,:));%染色体的适应度
end

%找到最好的适应度
[bestfitness bestindex]=min(fitness);
zbest=pop(bestindex,:);%全局最佳
gbest=pop;%个体最佳
fitnessgbest=fitness;%个体最佳适应度值
fitnesszbest=bestfitness;%全局最佳适应度值

%%进行迭代寻优
for i=1:maxgen
    i;
    for j=1:sizepop
        %进行速度更新
        V(j,:)= V(j,:)+c1*rand*(gbest(j,:)-pop(j,:))+c2*rand*(zbest-pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax;
        V(j,find(V(j,:)<Vmin))=Vmin;
        %进行种群更新
        pop(j,:)=pop(j,:)+0.5*V(j,:);
        pop(j,find(pop(j,:)>popmax))=popmax;
        pop(j,find(pop(j,:)<popmin))=popmin;
        
        %自适应性变异
        if rand>0.8
            k=ceil(par_num*rand);
            pop(j,k)=rand;
        end
        
        %适应度及边界限制条件
        if 0.072*pop(j,1)+0.063*pop(j,2)+0.057*pop(j,3)+0.05*pop(j,4)+0.032*pop(j,5)+...
            0.0442*pop(j,6)+0.0675*pop(j,7)<=264.4
             if 128*pop(j,1)+78.1*pop(j,2)+64.1*pop(j,3)+43*pop(j,4)+58.1*pop(j,5)+...
            36.9*pop(j,6)+50.5*pop(j,7)<=69719
            fitness(j)=fun(pop(j,:));
             end
        end
        %个体最优更新
        if fitness(j)<fitnessgbest(j)
            gbest(j,:)=pop(j,:);
            fitnessgbest(j)=fitness(j);
        end
        %群体最优更新
        if fitness(j)<fitnesszbest
            zbest=pop(j,:);
            fitnesszbest=fitness(j);
        end
    end
    yy(i)=fitnesszbest;
end

zbest
plot(yy);
title(['适应度曲线''终止代数=',num2str(maxgen)]);
xlabel('进化代数');
ylabel('适应度');
        
            
        
        
        
        



