function [par,par_best]=update_par(par,par_best)
    
    %Px=Px+Pv*t,����t=1,PxΪ��ǰ���ӵ�λ�ã�PvΪ��ǰ���ӵ��ٶ�
    par.x=par.x+par.vx;   
    par.y=par.x+par.vy;   
    
    par.fit=compute_fit(par);    %���㵱ǰ������Ӧ��
    
    %Pv=Pv+(c1*rand*(Gx-Px))+(c2*rand*(PBx-Px))
    %����c1,c2Ϊ��������
    %GxΪ���������Ӧ�����ӵ�λ��
    %PBxΪ��ǰ���ӵ����λ��
    c1=2;
    c2=2;
    par.vx=par.vx+c1*rand()*(par_best.x-par.x)+c2*rand()*(par.bestx-par.x);   
    par.vy=par.vy+c1*rand()*(par_best.y-par.y)+c2*rand()*(par.besty-par.y);
 
    if  par.fit>par.bestfit      %�����ǰ������Ӧ��Ҫ���ڵ�ǰ���������Ӧ��
        par.bestfit=par.fit;    %����µ�ǰ���������Ӧ��
        par.bestx=par.x;        %���µ�ǰ�������λ��
        par.besty=par.y;
        if par.bestfit>par_best.fit     %�����ǰ���������Ӧ�Ⱥ������������Ӧ��
            par_best.fit=par.bestfit;   %��������������Ӧ��
            par_best.x=par.x;           %�����������λ��
            par_best.y=par.y;
        end
    end

end