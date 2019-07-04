%连通图中各顶点间最短距离的计算
function D=shortdis(W)
%W权值矩阵（各点距离，不相邻为inf）
%对于W(i,j)，如果两顶点间存在弧，则为弧的初值，否则为inf；
%当i=j时，W(i,j)=0

%输入示例
%W=[0 1 3 4;1 0 2 inf;3 2 0 5;4 inf 5 0];
n=length(W);
D=W;
m=1;
while m<=n
    for i=1:n
        for j=1:n
            if D(i,j)>D(i,m)+D(m,j)
                D(i,j)=D(i,m)+D(m,j);   %距离进行更新
            end
        end
    end
    m=m+1;
end
D;
