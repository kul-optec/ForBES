function opt = Process_Options(opt)

% General options

if ~isfield(opt, 'tol') || isempty(opt.tol), opt.tol = 1e-8; end
if ~isfield(opt, 'term') || isempty(opt.term), opt.customTerm = false;
else opt.customTerm = true; end
if ~isfield(opt, 'record') || isempty(opt.record), opt.toRecord = false;
else opt.toRecord = true; end
if ~isfield(opt, 'maxit') || isempty(opt.maxit), opt.maxit = 10000; end
if ~isfield(opt, 'beta') || isempty(opt.beta), opt.beta = 0.05; end
if ~isfield(opt, 'variant'), opt.variant = ''; end
if ~isfield(opt, 'display') || isempty(opt.display), opt.display = 2; end
if ~isfield(opt, 'report') || isempty(opt.report), opt.report = 1; end
if ~isfield(opt, 'useHessian') || isempty(opt.useHessian), opt.useHessian = false; end
if ~isfield(opt, 'metric') || isempty(opt.metric), opt.metric = @(x) x; end

% Methods (directions) options

if ~isfield(opt, 'modBroyden') || isempty(opt.modBroyden), opt.modBroyden = 3; end
if ~isfield(opt, 'deltaCurvature') || isempty(opt.deltaCurvature), opt.deltaCurvature = 1e-6; end
if ~isfield(opt, 'thetaBar') || isempty(opt.thetaBar), opt.thetaBar = 1e-4; end
if ~isfield(opt, 'initialScaling') || isempty(opt.initialScaling), opt.initialScaling = 0; end
if ~isfield(opt, 'memory') || isempty(opt.memory), opt.memory = 10; end

opt.optsL.UT = true;
opt.optsL.TRANSA = true;
opt.optsU.UT = true;

if strcmp(opt.variant, 'fast'), opt.fast = 1;
else opt.fast = 0; end

% Sets default solver if not specified

if ~isfield(opt, 'solver') || isempty(opt.solver)
    opt.solver = 'zerofpr';
end
opt.solverfun = str2func(opt.solver);

% Sets default method if not specified

solver2method = containers.Map({'fbs', 'classical', 'minfbe', 'zerofpr', 'amls'}, ...
                               {'',    'lbfgs',     'lbfgs',  'lbfgs',   'lbfgs'});
if ~isfield(opt, 'method') || isempty(opt.method)
   opt.method = solver2method(opt.solver);
end
opt.methodfun = str2func(strcat('Direction_', lower(opt.method)));
if ~isfield(opt, 'memopt'), opt.memopt = 1; end

% Sets default line-search if not specified

if strcmp(opt.solver, 'classical')
    method2linesearch = containers.Map( ...
        {'sd',     'bfgs',       'lbfgs',      'cg-desc',    'cg-prp',     'cg-dyhs',    'bb',                 'broyden',    'lbroyden',   'rbroyden'}, ...
        {'armijo', 'lemarechal', 'lemarechal', 'lemarechal', 'lemarechal', 'lemarechal', 'nonmonotone-armijo', 'lemarechal', 'lemarechal', 'lemarechal'});
elseif strcmp(opt.solver, 'minfbe')
    method2linesearch = @(s) 'backtracking';
elseif strcmp(opt.solver, 'zerofpr')
    method2linesearch = @(s) 'backtracking';
elseif strcmp(opt.solver, 'amls')
    method2linesearch = @(s) 'backtracking';
else
    method2linesearch = @(s) '';
end
if ~isfield(opt, 'linesearch') || isempty(opt.linesearch)
    opt.linesearch = method2linesearch(opt.method);
end

% Wrap up a string describing the algorithm

opt.name = opt.solver;
if ~isempty(opt.variant) > 0, opt.name = strcat([opt.name, ', ', opt.variant]); end
if ~isempty(opt.method) > 0, opt.name = strcat([opt.name, ', ', opt.method]); end
if ~isempty(opt.linesearch) > 0, opt.name = strcat([opt.name, ', ', opt.linesearch]); end
opt.processed = true;

end
