function [fstar,W,flag]=maxflow(C)
%CΪ����·�����������
%fΪ�����
%WΪ�������
%flagΪ���, �ɴ˿ɵ���С�����ŵ�Ϊһ�飬δ����ŵ�Ϊһ��

%����ʾ��
%C=[0,13,9,0,0,0,0;0,0,0,6,5,0,0;0,4,0,5,6,5,0;0,0,0,0,5,4,4;...
%    0,0,0,0,0,0,9;0,0,0,0,0,0,10;0,0,0,0,0,0,0];
n=size(C,1);
%ȡ��ʼ������fΪ����
for i=1:n
    for(j=1:n)
        f(i,j)=0;
    end;
end
%��ʼ������i�ı��(flag(i),delta(i))
for i=1:n
    flag(i)=0;
    delta(i)=0;
end 
while(1)
    %������vs��ţ���"n+1"����"-"
    %��ʱvs�Ǳ�Ŷ�δ���ĵ㣬�������Ϊδ��ŵ�
    flag(1)=n+1;
    delta(1)=Inf; 
    while(1)
        %��Ź���,label=1˵��������δ��ţ�label=0˵���������ѱ��
        label=1; 
        for i=1:n
            %ѡ��һ���ѱ�ŵĵ�vi���м��
            if(flag(i)) 
           %��һ��δ��Ŷ���vi�����ĵ�vj���м�飬��ȷ���Ƿ��vj���б��
                for j=1:n
                    if (C(i,j)~=0|C(j,i)~=0) %�ж�vj��vi�Ƿ����
                        %��δ��ŵĵ�vj����(vi,vj)Ϊ�Ǳ��ͻ�ʱ
                        if(flag(j)==0&f(i,j)<C(i,j)) 
                            flag(j)=i;
                            delta(j)=C(i,j)-f(i,j);
                            label=0;
                            delta(j)=min(delta(j),delta(i));
                        %��δ��ŵĵ�vj, ��(vj,vi)Ϊ��������ʱ
                        elseif(flag(j)==0&(C(i,j)~=0|C(j,i)~=0)&f(j,i)>0) 
                            flag(j)=-i;
                            delta(j)=f(j,i);
                            label=0;
                            delta(j)=min(delta(j),delta(i));
                        end;
                    end;
                end;
            end;
        end
        %����յ�vt�ѱ���Ż���̵��м�����޷����,����ֹ��Ź���
        if(flag(n)|label)
            break;
        end;
    end
     %�յ�vtδ�����, f��Ϊ�����, �㷨��ֹ
    if(label)
        break;
    end
    %ȷ����������deltatΪ������
    deltat=delta(n);
    t=n; 
    %�����������,
    while(1)
        if(flag(t)>0)
            f(flag(t),t)=f(flag(t),t)+deltat; %ǰ�򻡵���
        elseif(flag(t)<0)
            f(flag(t),t)=f(flag(t),t)-deltat; %���򻡵���
        end 
        %��t�ı��Ϊvsʱ, ��ֹ��������
        if(flag(t)==1)
            for i=1:n
                flag(i)=0;
                delta(i)=0; 
            end;
            break;
        end 
        %��������ǰһ�λ��ϵ���f
        t=flag(t);
    end;
end; 
W=0;
fstar=f;
%�����������
for j=1:n
    W=W+f(1,j);
end
