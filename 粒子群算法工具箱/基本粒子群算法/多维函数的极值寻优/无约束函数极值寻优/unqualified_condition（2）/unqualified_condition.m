clc;
clear;
close all;

%������ʼ��
c1=1.49445;
c2=1.49445;
maxg=200;%��������
sizepop=20;%��Ⱥ��ģ

%�����ٶȺ���Ⱥ���±߽�ֵ
Vmax=1;
Vmin=-1;
popmax=5;
popmin=-5;

%������ʼ���Ӽ����ٶ�
for i=1:sizepop
    %�������һ����Ⱥ
    pop(i,:)=5*rands(1,2);%��ʼ����Ⱥ
    V(i,:)=rands(1,2);%��ʼ���ٶ�
    fitness(i)=fun(pop(i,:));%����Ⱦɫ�����Ӧ��
end
[bestfitness bestindex]=min(fitness);%Ѱ����Сֵ
zbest=pop(bestindex,:);%ȫ������
gbest=pop;             %��������
fitnessgbest=fitness;  %���������Ӧ��ֵ
fitnesszbest=bestfitness;%ȫ�������Ӧ��ֵ

%%����Ѱ��
for i=1:maxg
    %������������������
    for j=1:sizepop
        %�ٶȸ���
        V(j,:)= V(j,:)+c1*rand*(gbest(j,:)-pop(j,:))+c2*rand*(zbest-pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax;
        V(j,find(V(j,:)<Vmin))=Vmin;
        
        %��Ⱥ����
        pop(j,:)=pop(j,:)+0.5*V(j,:);
        pop(j,find(V(j,:)>popmax))=popmax;
        pop(j,find(V(j,:)<popmin))=popmin;
        
        %����Ӧ����
        if rand>0.8
            k=ceil(2*rand);
            pop(j,k)=rand;
        end
        
        %��Ӧ��ֵ
        fitness(j)=fun(pop(j,:));
        
        %�����������
        if fitness(j)<fitnessgbest(j)
            gbest(j,:)=pop(j,:);
            fitnessgbest(j)=fitness(j);
        end
        
        %Ⱥ���������
        if fitness(j)<fitnesszbest
            zbest=pop(j,:);
            fitnesszbest=fitness(j);
        end
    end
    yy(i)=fitnesszbest;
end
%%�������
plot(yy,'Linewidth',2);
title(['��Ӧ������''��ֹ����=',num2str(maxg)]);
grid on
xlabel('��������');
ylabel('��Ӧ��');

zbest   %��Ѹ���ֵ

    
