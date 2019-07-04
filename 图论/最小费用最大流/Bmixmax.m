function [f,W,zwf]=Bmixmax(b,C)
%b费用，C流量
%f为最小费用最大流
%W为最小费用最大流量
%zwf为最小费用

%输入示例
%C=[0 10 8 0 0;0 0 0 2 7;0 5 0 10 0;0 0 0 0 4;0 0 0 0 0];
%b=[0 4 1 0 0;0 0 0 6 1;0 2 0 3 0;0 0 0 0 2;0 0 0 0 0];
n=size(C,1);
W=0;
pre_bound=Inf; 
%取初始可行流f为零流
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
    %构造有向赋权图
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
    %用Ford 算法求最短路, 赋初值
    for i=2:n
        p(i)=Inf;
        s(i)=i;
    end     
    %求有向赋权图中vs到vt的最短路
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
    end %求最短路的Ford算法结束    
    if(p(n)==Inf)
        break;
    end %不存在vs到vt的最短路, 算法终止       
    deltat=Inf;     %deltat表示调整量
    t=n;     
    while(1) %计算调整量
        if(a(s(t),t)>0)
            dvtt=C(s(t),t)-f(s(t),t); %前向弧调整量
        elseif(a(s(t),t)<0)
            dvtt=f(t,s(t));end %后向弧调整量
        if(deltat>dvtt)deltat=dvtt;
        end
        if(s(t)==1)
            break;
        end     %当t的标号为vs时, 终止计算调整量
        t=s(t);
    end    %调整前一段弧上的流f
    pd=0;
    %如果最大流量大于或等于预定的流量值
    if(W+deltat>=pre_bound)
        deltat=pre_bound-W;
        pd=1;
    end
    t=n;   
    %调整过程
    while(1) 
        if(a(s(t),t)>0)
            f(s(t),t)=f(s(t),t)+deltat; %前向弧调整
        elseif(a(s(t),t)<0)
            f(t,s(t))=f(t,s(t))-deltat; %后向弧调整
        end         
        %当t的标号为vs时, 终止调整过程
        if(s(t)==1)
            break;
        end 
        t=s(t);
    end    
    %如果最大流量达到预定的流量值
    if(pd)
        break;
    end    
    %计算最大流量
    W=0; 
    for j=1:n
        W=W+f(1,j);
    end;
end 
%计算最小费用
zwf=0;
for(i=1:n)
    for(j=1:n)
        zwf=zwf+b(i,j)*f(i,j);
    end;
end
