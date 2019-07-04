%MenteCarlo求区域面积
clear
N=10000;
n=100;
for j=1:n
    k=0;
    for i=1:N
        a=12*rand(1,2)-6;%[-6,6]
        x(i)=a(1,1);
        y(i)=a(1,2);
        a1=(x(i)^2)/9+(y(i)^2)/36;
        a2=(x(i)^2)/36+y(i)^2;
        a3=(x(i)-2)^2+(y(i)+1)^2;
        if a1<1
            if a2<1
                if a3<9
                    k=k+1;
                end
            end
        end
    end
    m(j)=12^2*k/N;
end
mj=mean(m)
