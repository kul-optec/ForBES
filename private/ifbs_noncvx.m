function out = ifbs_noncvx(prob, opt, varargin)

MAXIMUM_Lf = 1e14;

% initialize output stuff

if opt.report
    residual = zeros(1, opt.maxit);
    objective = zeros(1, opt.maxit);
    ts = zeros(1, opt.maxit);
    % initialize operations counter
    ops = FBOperations();
else
    ops = [];
end

msgTerm = '';
record = [];

% display stuff

if opt.display >= 2
    fprintf('\n%s', opt.name);
    fprintf('\n%6s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.');
end

% get Lipschitz constant & adaptiveness

[Lf, adaptive] = prob.Get_Lipschitz(opt);

if adaptive
    warning('Lipschitz constant is not accurate')
end

% set stepsize, initialize vectors

gam = (0.99999-opt.beta)/Lf;

cache_x0 = FBCache(prob, prob.x0, gam, ops);
cache_x1 = FBCache(prob, prob.x0, gam, ops);
cache_0 = cache_x0;

tic0 = tic();

g_x1 = inf;

for it = 1:opt.maxit
    
    x0 = cache_x0.Get_Point();
    x1 = cache_x1.Get_Point();
    
    y1 = cache_x1.Get_GradStep();
    f_x1 = cache_x1.Get_f();
    w1 = y1 + (opt.beta/2)*(x1-x0);
    
    cache_w1 = FBCache(prob, w1, gam, ops);
    
    [x1, g_x1] = cache_w1.Get_ProxStep(w1);
    
    cache_x0 = cache_x1;
    cache_x1 = FBCache(prob, x1, gam, ops);
    
    % record values, stopping criterion
    
    if opt.toRecord
        record = [record, opt.record(prob, it, cache_0, cache_x1)];
    end
    
    res_curr = norm(x0-x1, 'inf')/cache_x1.Get_Gamma();
    
    if opt.report
        if it == 1
            objective(1, it) = inf;
        else
            objective(1, it) = f_x1 + g_x1;
        end
        residual(1, it) = res_curr;
        ts(1, it) = toc(tic0);
    end
    
    if res_curr <= opt.tol
        msgTerm = 'reached optimum (up to tolerance)';
        flagTerm = 0;
        break;
    end
    
    %%%
    
    % display stuff

    if opt.display == 1
        Util_PrintProgress(it);
    elseif (opt.display == 2 && mod(it,100) == 0) || opt.display >= 3
        obj_curr = f_x1 + g_x1;
        fprintf('%6d %7.4e %7.4e %7.4e\n', it, cache_x1.Get_Gamma(), res_curr, obj_curr);
    end

end

time = toc(tic0);

if it == opt.maxit
    msgTerm = [msgTerm, 'exceeded maximum iterations'];
    flagTerm = 1;
end

if opt.display == 1
    Util_PrintProgress(it, flagTerm);
elseif opt.display >= 2
    obj_curr = f_x1 + g_x1;
    fprintf('%6d %7.4e %7.4e %7.4e\n', it, cache_x1.Get_Gamma(), res_curr, obj_curr);
end

% pack up results

out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
out.x = x1;
out.iterations = it;
out.operations = ops;
if opt.report
    out.residual = residual(1, 1:it);
    out.objective = objective(1, 1:it);
    out.ts = ts(1, 1:it);
end
out.record = record;
out.gam = cache_x1.Get_Gamma();
out.adaptive = adaptive;
out.time = time;
