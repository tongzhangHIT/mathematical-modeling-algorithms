%��ͨͼ�и��������̾���ļ���
function D=shortdis(W)
%WȨֵ���󣨸�����룬������Ϊinf��
%����W(i,j)��������������ڻ�����Ϊ���ĳ�ֵ������Ϊinf��
%��i=jʱ��W(i,j)=0

%����ʾ��
%W=[0 1 3 4;1 0 2 inf;3 2 0 5;4 inf 5 0];
n=length(W);
D=W;
m=1;
while m<=n
    for i=1:n
        for j=1:n
            if D(i,j)>D(i,m)+D(m,j)
                D(i,j)=D(i,m)+D(m,j);   %������и���
            end
        end
    end
    m=m+1;
end
D;
