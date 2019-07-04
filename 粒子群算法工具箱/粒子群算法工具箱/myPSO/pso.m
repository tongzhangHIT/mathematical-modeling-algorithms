function [OUT,varargout]=pso(functname,D,varargin)
%% pso.m is changed from pso_Trelea_vectorized.m by C.X
%% 10/27/2009
%%
%% Brian Birge
%% Rev 3.3
%% 2/18/06
%%
%% pso.m : a generic particle swarm optimizer
%% to find the minimum or maximum of any MISO matlab function
%%
%% Usage:
%%   [optOUT]=pso(functname,D)
%% or:
%%   [optOUT,tr,bestpos]=pso(functname,D,mv,VarRange,minmax,PSOparams,PSOseedValue)
%%
%% Inputs:
%%    functname - �Ż���������������M�ļ���д
%%        D     - ���������������Ż�������ά��
%%
%% Optional Inputs:
%%       mv     - ��������ٶȣ�һ��ȡ�仯��Χ��10%~20%
%%    VarRange  - �����仯��Χ����ɾ��󣩣�����һ����ʽΪ��
%%                [ min1 max1
%%                  min2 max2
%%                     ...
%%                  minD maxD ]
%%     minmax   - �Ż�Ŀ�꣬��Ҫ������ֵ������Сֵ:
%%                = 0 ��ȡ��Сֵ��Ĭ�ϣ�
%%                = 1 ��ȡ���ֵ
%%                = 2 ��P(12)��Ӧ
%%    PSOparams - PSO����:
%%                P(1)-��ʾ����������
%%                P(2)-������������ʼ�����ٸ�����
%%                P(3)-����λ�õĴ�����
%%                     Ĭ��Ϊ3��������Խ�޺������ٶȷ���
%%                P(4)-�㷨�ļ��ٶȲ�����Ӱ��ֲ�����ֵ����˵2�Դ����
%%                     �����˵����ͦ�õ�ѡ������Ĭ������Ϊ2
%%                P(5)-�㷨�ļ��ٶȲ�����Ӱ��ȫ������ֵ����˵2�Դ����
%%                     �����˵����ͦ�õ�ѡ������Ĭ������Ϊ2
%%                P(6)-��ʼʱ�̵ļ�Ȩֵ������������ٶȺ��������ȣ�Ĭ��
%%                     ����Ϊ0.9
%%                P(7)-����ʱ�̵ļ�Ȩֵ������������ٶȺ��������ȣ�Ĭ��
%%                     ����Ϊ0.4
%%                P(8)-ָ���ĵ���������������ֵʱ����Ȩֵȡ����С(�������0.4)
%%                P(9)-������ֹ�㷨����ֵ�������������ε����ж�Ӧ����Ⱥ
%%                     ����ֵС�ڴ���ֵʱ���㷨ֹͣ
%%                P(10)-������ֹ�㷨����ֵ��������P(10)�ε����к�������
%%                      ��ֵ��Ȼû�б仯�����Ƴ�����
%%                P(11)-����˵���Ż��������ȡNaNʱ��ʾΪ��Լ���µ��Ż�
%%                      ����(��û�и���Լ������)
%%                P(12)-����ָ������PSO�����ͣ�0��ʾͨ����PSO�㷨
%%                P(13)-����˵���Ƿ����û����в������ӣ�0��ʾ���������
%%                      �ӣ�1��ʾ�û����в�������
%% PSOseedValue - ���û����в��������ӣ���P(13)=0,��������ò���
%%
%% Outputs:
%%      optOUT  - �õ�������Ż����������Ż�����ֵ����һ����ʽΪ��
%%                [ bestin1
%%                  bestin2
%%                    ...
%%                  bestinD
%%                  bestOUT ]
%%
%% Optional Outputs:
%%        tr    - ÿһ�����ӵ�ȫ������ֵ���������ӵ����й켣
%%     bestpos  - ÿһ�����ӵ�ȫ�����Ų�����ȫ������ֵ����һ����ʽΪ��
%%                [ gbest1 gbestval1
%%                  gbest2 gbestval2
%%                       ...
%%                  gbestP(1) gbestvalP(1) ]


tic;
Handle = waitbar(0,'Please wait...');
pause(0.1);
rand('state',sum(100*clock));

if nargin < 2
   error('Not enough arguments.');
end

%---------------------------PSO PARAMETERS
if nargin == 2      % only specified functname and D
   VRmin=ones(D,1)*-100; 
   VRmax=ones(D,1)*100;    
   VR=[VRmin,VRmax];
   minmax = 0;
   P = [];
   mv = 4;  
elseif nargin == 3  % specified functname, D, and mv
   VRmin=ones(D,1)*-100; 
   VRmax=ones(D,1)*100;    
   VR=[VRmin,VRmax];
   minmax = 0;
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end
   P = [];  
elseif nargin == 4  % specified functname, D, mv, Varrange
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end
   VR=varargin{2}; 
   minmax = 0;
   P = [];   
elseif nargin == 5  % Functname, D, mv, Varrange, and minmax
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end    
   VR=varargin{2};
   minmax=varargin{3};
   P = [];
elseif nargin == 6  % Functname, D, mv, Varrange, minmax, and psoparams
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end    
   VR=varargin{2};
   minmax=varargin{3};
   P = varargin{4}; % psoparams
elseif nargin == 7  % Functname, D, mv, Varrange, minmax, and psoparams, PSOseedValue
   mv=varargin{1};
   if isnan(mv)
       mv=4;
   end    
   VR=varargin{2};
   minmax=varargin{3};
   P = varargin{4}; % psoparams 
   PSOseedValue = varargin{5};
else    
   error('Wrong # of input arguments.');
end

% sets up default pso params
Pdef = [2000 24 3 2 2 0.9 0.4 1500 1e-6 250 NaN 0 0];
Plen = length(P);
P    = [P,Pdef(Plen+1:end)];
me          = P(1);
ps          = P(2);
posmaskmeth = P(3);
ac1         = P(4);
ac2         = P(5);
iw1         = P(6);
iw2         = P(7);
iwe         = P(8);
ergrd       = P(9);
ergrdep     = P(10);
errgoal     = P(11);
trelea      = P(12);
PSOseed     = P(13);
%---------------------------PSO PARAMETERS setting ends

%---------------------------error checking
if ((minmax==2) & isnan(errgoal))
    error('minmax= 2, errgoal= NaN: choose an error goal or set minmax to 0 or 1');
end
if ( (PSOseed==1) & ~exist('PSOseedValue') )
    error('PSOseed flag set but no PSOseedValue was input');
end
if exist('PSOseedValue')
    tmpsz=size(PSOseedValue);
    if D < tmpsz(2)
        error('PSOseedValue column size must be D or less');
    end
    if ps < tmpsz(1)
        error('PSOseedValue row length must be # of particles or less');
    end
end
%---------------------------error checking ends

% take care of setting max velocity and position params here
if length(mv)==1
    velmaskmin = -mv*ones(ps,D);     % min vel
    velmaskmax =  mv*ones(ps,D);     % max vel
elseif length(mv)==D     
    velmaskmin = repmat(forcerow(-mv),ps,1); % min vel
    velmaskmax = repmat(forcerow( mv),ps,1); % max vel
else
    error('Max vel must be either a scalar or same length as prob dimension D');
end
posmaskmin  = repmat(VR(1:D,1)',ps,1);  % min pos
posmaskmax  = repmat(VR(1:D,2)',ps,1);  % max pos

%---------------------------initialize
pos(1:ps,1:D) = normmat(rand([ps,D]),VR',1);
if PSOseed == 1      % initial positions user input, see comments above
   tmpsz                      = size(PSOseedValue);
   pos(1:tmpsz(1),1:tmpsz(2)) = PSOseedValue;  
end
vel(1:ps,1:D) = normmat(rand([ps,D]),[forcecol(-mv),forcecol(mv)]',1);

% initial pbest positions vals
pbest   = pos;
out     = feval(functname,pos);
pbestval=out;   % initially, pbest is same as pos

% assign initial gbest here also (gbest and gbestval)
if minmax==1
   % this picks gbestval when we want to maximize the function
   [gbestval,idx1] = max(pbestval);
elseif minmax==0
   % this works for straight minimization
   [gbestval,idx1] = min(pbestval);
elseif minmax==2
   % this works when you know target but not direction you need to go
   % good for a cost function that returns distance to target that can be either
   % negative or positive (direction info)
   [temp,idx1] = min((pbestval-ones(size(pbestval))*errgoal).^2);
   gbestval    = pbestval(idx1);
end
gbest      = pbest(idx1,:);  % this is gbest position

% preallocate variables for speed up
tr = ones(1,me)*NaN;
% preallocate a variable to keep track of gbest for all iters
bestpos    = zeros(me,D+1)*NaN;
%---------------------------initialize ends

cut = 0;
iwt(1) = iw1;
generations = 1;

%---------------------------start epoch loop (iterations)
for i=1:me  % start epoch loop (iterations)

% find particles where we have new pbest, depending on minmax choice 
% then find gbest and gbestval
    if minmax == 0
        [tempi]            = find(pbestval>=out); % new min pbestvals
        pbestval(tempi,1)  = out(tempi);   % update pbestvals
        pbest(tempi,:)     = pos(tempi,:); % update pbest positions
        [iterbestval,idx1] = min(pbestval);
        if gbestval >= iterbestval
            gbestval = iterbestval;
            gbest    = pbest(idx1,:);
        end
    elseif minmax == 1
        [tempi,dum]        = find(pbestval<=out); % new max pbestvals
        pbestval(tempi,1)  = out(tempi,1); % update pbestvals
        pbest(tempi,:)     = pos(tempi,:); % update pbest positions
        [iterbestval,idx1] = max(pbestval);
        if gbestval <= iterbestval
            gbestval = iterbestval;
            gbest    = pbest(idx1,:);
        end
    elseif minmax == 2  % this won't work as it is, fix it later
        egones            = errgoal*ones(ps,1); % vector of errgoals
        sqrerr2           = ((pbestval-egones).^2);
        sqrerr1           = ((out-egones).^2);
        [tempi,dum]       = find(sqerr1 <= sqrerr2); % find particles closest to targ
        pbestval(tempi,1) = out(tempi,1); % update pbestvals
        pbest(tempi,:)    = pos(tempi,:); % update pbest positions
        sqrerr            = ((pbestval-egones).^2); % need to do this to reflect new pbests
        [temp,idx1]       = min(sqrerr);
        iterbestval       = pbestval(idx1);
        if (iterbestval-errgoal)^2 <= (gbestval-errgoal)^2
            gbestval = iterbestval;
            gbest    = pbest(idx1,:);
        end
    end
    
% PSOPSOPSOPSOPSOPSO-->
    % get new velocities, positions (this is the heart of the PSO algorithm)     
    % each epoch get new set of random numbers
    rannum1 = rand([ps,D]); % for Trelea and Clerc types
    rannum2 = rand([ps,D]);       
    if trelea == 2    
    % from Trelea's paper, parameter set 2
        vel = 0.729.*vel...                              % prev vel
              +1.494.*rannum1.*(pbest-pos)...            % independent
              +1.494.*rannum2.*(repmat(gbest,ps,1)-pos); % social  
    elseif trelea == 1
    % from Trelea's paper, parameter set 1                     
        vel = 0.600.*vel...                              % prev vel
              +1.700.*rannum1.*(pbest-pos)...            % independent
              +1.700.*rannum2.*(repmat(gbest,ps,1)-pos); % social 
    elseif trelea ==3
    % Clerc's Type 1" PSO
        vel = chi*(vel...                                % prev vel
              +ac1.*rannum1.*(pbest-pos)...              % independent
              +ac2.*rannum2.*(repmat(gbest,ps,1)-pos)) ; % social          
    else
    % common PSO algo with inertia wt 
    % get inertia weight, just a linear funct w.r.t. epoch parameter iwe
        if i<=iwe
            iwt(i) = ((iw2-iw1)/(iwe-1))*(i-1)+iw1;
        else
            iwt(i) = iw2;
        end
        % random number including acceleration constants
        ac11 = rannum1.*ac1;    % for common PSO w/inertia
        ac22 = rannum2.*ac2;
        vel = iwt(i).*vel...                             % prev vel
              +ac11.*(pbest-pos)...                      % independent
              +ac22.*(repmat(gbest,ps,1)-pos);           % social                  
    end
    
    % limit velocities here using masking
    vel = ( (vel <= velmaskmin).*velmaskmin ) + ( (vel > velmaskmin).*vel );
    vel = ( (vel >= velmaskmax).*velmaskmax ) + ( (vel < velmaskmax).*vel );      
    % update new position (PSO algo)    
    pos = pos + vel;
    % position masking, limits positions to desired search space
    % method: 0) no position limiting, 1) saturation at limit,
    %         2) wraparound at limit , 3) bounce off limit
    minposmask_throwaway = pos <= posmaskmin;  % these are psXD matrices
    minposmask_keep      = pos >  posmaskmin;     
    maxposmask_throwaway = pos >= posmaskmax;
    maxposmask_keep      = pos <  posmaskmax;
     
    if     posmaskmeth == 1
    % this is the saturation method
        pos = ( minposmask_throwaway.*posmaskmin ) + ( minposmask_keep.*pos );
        pos = ( maxposmask_throwaway.*posmaskmax ) + ( maxposmask_keep.*pos );      
    elseif posmaskmeth == 2
    % this is the wraparound method
        pos = ( minposmask_throwaway.*posmaskmax ) + ( minposmask_keep.*pos );
        pos = ( maxposmask_throwaway.*posmaskmin ) + ( maxposmask_keep.*pos );                
    elseif posmaskmeth == 3
    % this is the bounce method, particles bounce off the boundaries with -vel      
        pos = ( minposmask_throwaway.*posmaskmin ) + ( minposmask_keep.*pos );
        pos = ( maxposmask_throwaway.*posmaskmax ) + ( maxposmask_keep.*pos );
        vel = (vel.*minposmask_keep) + (-vel.*minposmask_throwaway);
        vel = (vel.*maxposmask_keep) + (-vel.*maxposmask_throwaway);
    end
% PSOPSOPSOPSOPSOPSO<--

    tr(i+1)          = gbestval;
    bestpos(i,1:D+1) = [gbest,gbestval];
    out        = feval(functname,pos);
    
    mmm=round(generations*100/me);
    waitbar(mmm/100,Handle,[num2str(mmm),'%finished...'],Handle);
    pause(0.1);
    generations = generations+1;
    
% check for stopping criterion based on speed of convergence to desired    
    tmp1 = abs(tr(i) - gbestval);
    if tmp1 > ergrd
        cut = 0;
    elseif tmp1 <= ergrd
        cut = cut+1;
        if cut >= ergrdep       
            break
        end
    end
    
% this stops if using constrained optimization and goal is reached
    if ~isnan(errgoal)
        if ((gbestval<=errgoal) & (minmax==0)) | ((gbestval>=errgoal) & (minmax==1))
            break
        end
    % this is stopping criterion for constrained from both sides    
        if minmax == 2
            if ((tr(i)<errgoal) & (gbestval>=errgoal)) | ((tr(i)>errgoal) ...
                    & (gbestval <= errgoal))        
                break              
            end
        end
    end

end
%---------------------------end epoch loop

% output & return
OUT=[gbest';gbestval];
varargout{1}=[tr(find(~isnan(tr)))];
varargout{2}=bestpos(1:i,1:D+1);

waitbar(100/100,Handle,[num2str(100),'%finished...'],Handle);
pause(0.1);
close(Handle);
toc;

return