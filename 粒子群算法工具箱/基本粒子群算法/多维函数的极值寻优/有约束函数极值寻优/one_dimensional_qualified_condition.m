clc
clear 
close all
warning off

%%�����ĳ�ʼ��
c1=1.49445;%---------�ɽ��и���-------
c2=1.49445;%---------�ɽ��и���-------
maxgen=200; %��������%---------�ɽ��и���-------
sizepop=200;%��Ⱥ��ģ%---------�ɽ��и���-------

Vmax=1;     %���Ӹ����ٶ�
Vmin=-1;

popmax=50;%��Ⱥ
popmin=-50;

par_num=7;%ά��

%%������ʼ���Ӻ��ٶ�
for i=1:sizepop
    pop(i,:)=1.*rands(1,par_num);%��ʼ����Ⱥ
    V(i,:)=1.*rands(1,par_num);%��ʼ���ٶ�
    %������Ӧ��
    fitness(i)=fun(pop(i,:));%Ⱦɫ�����Ӧ��
end

%�ҵ���õ���Ӧ��
[bestfitness bestindex]=min(fitness);
zbest=pop(bestindex,:);%ȫ�����
gbest=pop;%�������
fitnessgbest=fitness;%���������Ӧ��ֵ
fitnesszbest=bestfitness;%ȫ�������Ӧ��ֵ

%%���е���Ѱ��
for i=1:maxgen
    i;
    for j=1:sizepop
        %�����ٶȸ���
        V(j,:)= V(j,:)+c1*rand*(gbest(j,:)-pop(j,:))+c2*rand*(zbest-pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax;
        V(j,find(V(j,:)<Vmin))=Vmin;
        %������Ⱥ����
        pop(j,:)=pop(j,:)+0.5*V(j,:);
        pop(j,find(pop(j,:)>popmax))=popmax;
        pop(j,find(pop(j,:)<popmin))=popmin;
        
        %����Ӧ�Ա���
        if rand>0.8
            k=ceil(par_num*rand);
            pop(j,k)=rand;
        end
        
        %��Ӧ�ȼ��߽���������
        if 0.072*pop(j,1)+0.063*pop(j,2)+0.057*pop(j,3)+0.05*pop(j,4)+0.032*pop(j,5)+...
            0.0442*pop(j,6)+0.0675*pop(j,7)<=264.4
             if 128*pop(j,1)+78.1*pop(j,2)+64.1*pop(j,3)+43*pop(j,4)+58.1*pop(j,5)+...
            36.9*pop(j,6)+50.5*pop(j,7)<=69719
            fitness(j)=fun(pop(j,:));
             end
        end
        %�������Ÿ���
        if fitness(j)<fitnessgbest(j)
            gbest(j,:)=pop(j,:);
            fitnessgbest(j)=fitness(j);
        end
        %Ⱥ�����Ÿ���
        if fitness(j)<fitnesszbest
            zbest=pop(j,:);
            fitnesszbest=fitness(j);
        end
    end
    yy(i)=fitnesszbest;
end

zbest
plot(yy);
title(['��Ӧ������''��ֹ����=',num2str(maxgen)]);
xlabel('��������');
ylabel('��Ӧ��');
        
            
        
        
        
        



