function re=compute_fit(par)
    x=par.x;
    y=par.y;
    if x<-0.5 || x>0.5 || y<-0.5 || y>0.5    %-----------���滻
        re=0;        %������Χ��Ӧ��Ϊ0
    else            %������Ӧ�Ȱ�Ŀ�꺯�����
        re=-20*exp(-0.2*sqrt((x^2+y^2)/2))-exp((cos(2*pi*x)+cos(2*pi*y))/2)+20+2.71928; 
    end
end