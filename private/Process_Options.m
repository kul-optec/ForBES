function [opt, lsopt] = Process_Options(opt)

% fill in missing options with defaults

if ~isfield(opt, 'tol') || isempty(opt.tol), opt.tol = 1e-8; end
if ~isfield(opt, 'term') || isempty(opt.term), opt.customTerm = false;
else opt.customTerm = true; end
if ~isfield(opt, 'record') || isempty(opt.record), opt.toRecord = false;
else opt.toRecord = true; end
if ~isfield(opt, 'maxit') || isempty(opt.maxit), opt.maxit = 10000; end
if ~isfield(opt, 'beta') || isempty(opt.beta), opt.beta = 0.05; end
if ~isfield(opt, 'Lf') || isempty(opt.Lf), opt.adaptive = 1; end
if ~isfield(opt, 'adaptive') || isempty(opt.adaptive), opt.adaptive = 0; end
if ~isfield(opt, 'variant'), opt.variant = ''; end
if ~isfield(opt, 'display') || isempty(opt.display), opt.display = 2; end
if ~isfield(opt, 'useHessian') || isempty(opt.useHessian), opt.useHessian = 0; end
if ~isfield(opt, 'bopt') || isempty(opt.bopt), opt.bopt = 0; end
if ~isfield(opt, 'memory') || isempty(opt.memory), opt.memory = 10; end
opt.optsL.UT = true;
opt.optsL.TRANSA = true;
opt.optsU.UT = true;

if strcmp(opt.variant, 'fast'), opt.fast = 1;
else opt.fast = 0; end

% sets default solver if not specified

if ~isfield(opt, 'solver') || isempty(opt.solver)
    opt.solver = 'zerofpr';
end
opt.solverfun = str2func(opt.solver);

% sets default method if not specified

solver2method = containers.Map({'fbs', 'classical', 'minfbe', 'zerofpr'}, ...
                               {'',    'lbfgs',     'lbfgs',  'lbfgs'});
if ~isfield(opt, 'method') || isempty(opt.method)
   opt.method = solver2method(opt.solver);
end
opt.methodfun = str2func(strcat('Direction_', lower(opt.method)));
if ~isfield(opt, 'memopt'), opt.memopt = 1; end

% sets default line-search if not specified

if strcmp(opt.solver, 'classical')
    method2linesearch = containers.Map( ...
        {'sd',     'bfgs',       'lbfgs',      'cg-desc',    'cg-prp',     'cg-dyhs',    'bb',                 'broyden',    'lbroyden',   'rbroyden'}, ...
        {'armijo', 'lemarechal', 'lemarechal', 'lemarechal', 'lemarechal', 'lemarechal', 'nonmonotone-armijo', 'lemarechal', 'lemarechal', 'lemarechal'});
elseif strcmp(opt.solver, 'minfbe')
    method2linesearch = @(s) 'backtracking';
elseif strcmp(opt.solver, 'zerofpr')
    method2linesearch = @(s) 'backtracking';
elseif strcmp(opt.solver, 'zerofpr_new')
    method2linesearch = @(s) 'backtracking';
else
    method2linesearch = @(s) '';
end
if ~isfield(opt, 'linesearch') || isempty(opt.linesearch)
    opt.linesearch = method2linesearch(opt.method);
end

% wrap up a string describing the algorithm

opt.name = opt.solver;
if length(opt.variant) > 0, opt.name = strcat([opt.name, ', ', opt.variant]); end
if length(opt.method) > 0, opt.name = strcat([opt.name, ', ', opt.method]); end
if length(opt.linesearch) > 0, opt.name = strcat([opt.name, ', ', opt.linesearch]); end
opt.processed = true;

end
