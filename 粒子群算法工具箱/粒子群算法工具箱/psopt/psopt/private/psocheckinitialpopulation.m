function [state,options,Aineq,bineq,Aeq,beq] = ...
    psocheckinitialpopulation(state,Aineq,bineq,Aeq,beq,LB,UB,nonlcon,...
    options)
% Checks initial population with respect to linear constraints. Requires
% optimization toolbox.

if exist('linprog','file') ~= 2
    msg = sprintf('Could not find a required function in Optimization ') ;
    msg = sprintf('%s Toolbox. Ignoring linear constraints ',msg) ;
    msg = sprintf('%s for initial population distribution and',msg) ;
    warning('pso:linearconstraints:missingtoolbox',...
        '%s setting constraint behavior to ''soft''.',msg)
    options.ConstrBoundary = 'soft' ;
    return
end

state.LinprogOptions = optimset('Simplex','off',...
                'LargeScale','off',...
                'Display','off') ;
state.OutOfBounds = false(options.PopulationSize,1) ;

hw = waitbar(0,'Finding feasible initial positions...') ;
for i = 1:size(state.Population,1)
    if (~isempty(Aineq) && ...
            any(Aineq*state.Population(i,:)' - bineq > options.TolCon)) ...
            || (~isempty(Aeq) && ...
            any(abs(Aeq*state.Population(i,:)' - beq) > options.TolCon))
        % Reposition the ith particle if it is outside of linear constraint
        % boundaries
        [newpoint,unused,exitflag] = ...
            linprog([],Aineq,bineq,Aeq,beq,LB,UB,...
            state.Population(i,:),...
            state.LinprogOptions) ;
        clear unused
        if exitflag == -2
            error('Problem is infeasible due to constraints')
        else
            state.Population(i,:) = reshape(newpoint,1,[]) ;
        end % if exitflag
    end % if any
    if ~isempty(nonlcon)
        [c,ceq] = nonlcon(state.Population(i,:)) ;
        % Basically, keep trying random points within PopInitRange until we
        % find one that satisfies the nonlinear constraints.
        if ~isempty(ceq) && i == 1
            msg = 'Nonlinear equality constraints cannot be solved by' ;
            msg = sprintf('%s PSO yet. Using only the inequality ',msg) ;
            warning('pso:constraints:nonlinear:equality',...
                '%s constraints for now...',msg)
        end
        while any(c > options.TolCon)
            state.Population(i,:) = options.PopInitRange(1,:) + ...
                rand(1,size(options.PopInitRange,2)) .* ...
                (options.PopInitRange(2,:) - options.PopInitRange(1,:)) ;
            c = nonlcon(state.Population(i,:)) ;
        end
    end % if ~isempty
    waitbar(i/options.PopulationSize,hw)
end % for i
if ishandle(hw), close(hw), end