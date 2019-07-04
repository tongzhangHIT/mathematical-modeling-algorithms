function [f,W,zwf]=Bmixmax(b,C)
%b���ã�C����
%fΪ��С���������
%WΪ��С�����������
%zwfΪ��С����

%����ʾ��
%C=[0 10 8 0 0;0 0 0 2 7;0 5 0 10 0;0 0 0 0 4;0 0 0 0 0];
%b=[0 4 1 0 0;0 0 0 6 1;0 2 0 3 0;0 0 0 0 2;0 0 0 0 0];
n=size(C,1);
W=0;
pre_bound=Inf; 
%ȡ��ʼ������fΪ����
for i=1:n
    for j=1:n
        f(i,j)=0;
    end;
end 
while(1)
    for i=1:n
        for j=1:n
            if j~=i
                a(i,j)=Inf;
            end
        end
    end    
    %��������Ȩͼ
    for i=1:n
        for j=1:n
            if(C(i,j)>0&f(i,j)==0)
                a(i,j)=b(i,j);
            elseif(C(i,j)>0&f(i,j)==C(i,j))
                a(j,i)=-b(i,j);
            elseif(C(i,j)>0)a(i,j)=b(i,j);
                a(j,i)=-b(i,j);
            end
        end
    end    
    %��Ford �㷨�����·, ����ֵ
    for i=2:n
        p(i)=Inf;
        s(i)=i;
    end     
    %������Ȩͼ��vs��vt�����·
    for k=1:n
        pd=1; 
        for i=2:n
            for j=1:n
                if(p(i)>p(j)+a(j,i))
                    p(i)=p(j)+a(j,i);
                    s(i)=j;
                    pd=0;
                end
            end
        end
        if(pd)
            break;
        end
    end %�����·��Ford�㷨����    
    if(p(n)==Inf)
        break;
    end %������vs��vt�����·, �㷨��ֹ       
    deltat=Inf;     %deltat��ʾ������
    t=n;     
    while(1) %���������
        if(a(s(t),t)>0)
            dvtt=C(s(t),t)-f(s(t),t); %ǰ�򻡵�����
        elseif(a(s(t),t)<0)
            dvtt=f(t,s(t));end %���򻡵�����
        if(deltat>dvtt)deltat=dvtt;
        end
        if(s(t)==1)
            break;
        end     %��t�ı��Ϊvsʱ, ��ֹ���������
        t=s(t);
    end    %����ǰһ�λ��ϵ���f
    pd=0;
    %�������������ڻ����Ԥ��������ֵ
    if(W+deltat>=pre_bound)
        deltat=pre_bound-W;
        pd=1;
    end
    t=n;   
    %��������
    while(1) 
        if(a(s(t),t)>0)
            f(s(t),t)=f(s(t),t)+deltat; %ǰ�򻡵���
        elseif(a(s(t),t)<0)
            f(t,s(t))=f(t,s(t))-deltat; %���򻡵���
        end         
        %��t�ı��Ϊvsʱ, ��ֹ��������
        if(s(t)==1)
            break;
        end 
        t=s(t);
    end    
    %�����������ﵽԤ��������ֵ
    if(pd)
        break;
    end    
    %�����������
    W=0; 
    for j=1:n
        W=W+f(1,j);
    end;
end 
%������С����
zwf=0;
for(i=1:n)
    for(j=1:n)
        zwf=zwf+b(i,j)*f(i,j);
    end;
end
