function [fstar,W,flag]=maxflow(C)
%C为各线路最大允许流量
%f为最大流
%W为最大流量
%flag为标号, 由此可得最小割，被标号的为一组，未被标号的为一组

%输入示例
%C=[0,13,9,0,0,0,0;0,0,0,6,5,0,0;0,4,0,5,6,5,0;0,0,0,0,5,4,4;...
%    0,0,0,0,0,0,9;0,0,0,0,0,0,10;0,0,0,0,0,0,0];
n=size(C,1);
%取初始可行流f为零流
for i=1:n
    for(j=1:n)
        f(i,j)=0;
    end;
end
%初始化顶点i的标号(flag(i),delta(i))
for i=1:n
    flag(i)=0;
    delta(i)=0;
end 
while(1)
    %给发点vs标号，用"n+1"代替"-"
    %这时vs是标号而未检查的点，其它点均为未标号点
    flag(1)=n+1;
    delta(1)=Inf; 
    while(1)
        %标号过程,label=1说明过程中未标号，label=0说明过程中已标号
        label=1; 
        for i=1:n
            %选择一个已标号的点vi进行检查
            if(flag(i)) 
           %对一切未标号而与vi关联的点vj进行检查，以确定是否对vj进行标号
                for j=1:n
                    if (C(i,j)~=0|C(j,i)~=0) %判断vj与vi是否关联
                        %对未标号的点vj，当(vi,vj)为非饱和弧时
                        if(flag(j)==0&f(i,j)<C(i,j)) 
                            flag(j)=i;
                            delta(j)=C(i,j)-f(i,j);
                            label=0;
                            delta(j)=min(delta(j),delta(i));
                        %对未标号的点vj, 当(vj,vi)为非零流弧时
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
        %如果收点vt已被标号或过程的中间点已无法标号,则终止标号过程
        if(flag(n)|label)
            break;
        end;
    end
     %收点vt未被标号, f已为最大流, 算法终止
    if(label)
        break;
    end
    %确定调整量，deltat为调整量
    deltat=delta(n);
    t=n; 
    %进入调整过程,
    while(1)
        if(flag(t)>0)
            f(flag(t),t)=f(flag(t),t)+deltat; %前向弧调整
        elseif(flag(t)<0)
            f(flag(t),t)=f(flag(t),t)-deltat; %后向弧调整
        end 
        %当t的标号为vs时, 终止调整过程
        if(flag(t)==1)
            for i=1:n
                flag(i)=0;
                delta(i)=0; 
            end;
            break;
        end 
        %继续调整前一段弧上的流f
        t=flag(t);
    end;
end; 
W=0;
fstar=f;
%计算最大流量
for j=1:n
    W=W+f(1,j);
end
