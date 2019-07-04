function [xOpt,fval] = psorunhybridfcn(fitnessfcn,xOpt,...
    Aineq,bineq,Aeq,beq,LB,UB,nonlcon,options)

fprintf('\nBest point before hybrid function: %s',...
    mat2str(xOpt,5))
fprintf('\n\nTurning over to hybrid function %s...\n\n',...
    func2str(options.HybridFcn))
hybridOptions = optimset('LargeScale','off') ;
% hybridOptions = optimset ;

if exist(func2str(options.HybridFcn),'file') ~= 2
    warning('pso:hybridfcn:nofile',...
        'Hybrid function %s cannot be found. Check toolboxes.',...
        func2str(options.HybridFcn))
    fval = fitnessfcn(xOpt) ;
    return
end

if strcmp(func2str(options.HybridFcn),func2str(@fmincon)) && ...
        all([isempty([Aineq,bineq]), isempty([Aeq,beq]), ...
        isempty([LB;UB]),isempty(nonlcon)])
    msg = sprintf('Warning: %s does not accept problems without any',...
        func2str(options.HybridFcn)) ;
    fprintf('%s constraints. Switching to fminunc.\n\n',msg)
    options.HybridFcn = @fminunc ;
elseif (strcmp(func2str(options.HybridFcn),func2str(@fminunc)) || ...
        strcmp(func2str(options.HybridFcn),func2str(@fminsearch))) && ...
        ~all([isempty([Aineq,bineq]), isempty([Aeq,beq]), ...
        isempty([LB;UB]),isempty(nonlcon)])
    msg = sprintf('Warning: %s does not accept problems with',...
        func2str(options.HybridFcn)) ;
    fprintf(0,'%s constraints. Switching to fmincon.\n\n',msg)
    options.HybridFcn = @fmincon ;
end

if strcmp(func2str(options.HybridFcn),func2str(@fmincon)) || ...
        strcmp(func2str(options.HybridFcn),func2str(@patternsearch))
    [xOpt,fval] = options.HybridFcn(fitnessfcn,xOpt,Aineq,bineq,...
        Aeq,beq,LB,UB,nonlcon,hybridOptions) ;
elseif strcmp(func2str(options.HybridFcn),func2str(@fminunc)) || ...
        strcmp(func2str(options.HybridFcn),func2str(@fminsearch))
    [xOpt,fval] = options.HybridFcn(fitnessfcn,xOpt,hybridOptions) ;
else
    warning('pso:hybridfcn:unrecognized',...
        'Unrecognized hybrid function. Ignoring for now.')
end