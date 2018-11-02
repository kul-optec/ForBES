function out = fbs(prob, opt, varargin)

MAXIMUM_Lf = 1e14;

% initialize output stuff

if opt.report
    residual = zeros(1, opt.maxit);
    objective = zeros(1, opt.maxit);
    ts = zeros(1, opt.maxit);
end

% initialize operations counter
ops = FBOperations();

msgTerm = '';
record = [];

% display stuff

if opt.display >= 2
    fprintf('\n%s', opt.name);
    fprintf('\n%6s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.');
end

% get Lipschitz constant & adaptiveness

[Lf, adaptive] = prob.Get_Lipschitz(opt);

% set stepsize, initialize vectors

gam = 1/Lf;
xk = prob.x0;
vk = prob.x0;

t0 = tic();

for it = 1:opt.maxit

    if opt.fast
        theta = 2/(it+1); % since it starts from 1
        yk = (1-theta)*xk+theta*vk;
    else
        yk = xk;
    end

    cache_yk = FBCache(prob, yk, gam, ops);

    if it == 1
        cache_0 = cache_yk;
    end

    hasGammaChanged = false;

    if adaptive
        [hasGammaChanged, ~] = cache_yk.Backtrack_Gamma(0.0);
        gam = cache_yk.Get_Gamma();
    end

    if opt.report
        objective(1, it) = cache_yk.Get_FBE();
        residual(1, it) = norm(cache_yk.Get_FPR(), 'inf')/cache_yk.Get_Gamma();
        ts(1, it) = toc(t0);
    end

    if opt.toRecord
        record = [record, opt.record(prob, it, cache_0, cache_yk)];
    end

    if ~hasGammaChanged
        if ~opt.customTerm
            if cache_yk.Check_StoppingCriterion(opt.tol)
                msgTerm = 'reached optimum (up to tolerance)';
                flagTerm = 0;
                break;
            end
        else
            %flagStop = opt.term(prob, it, cache_0, cache_yk);
            flagStop = opt.term(prob, it, cache_0, cache_yk, opt.hack);
            if (adaptive == 0 || it > 1) && flagStop
                msgTerm = 'reached optimum (custom criterion)';
                flagTerm = 0;
                break;
            end
        end
    end
    if prob.Lf >= MAXIMUM_Lf
        msgTerm = 'L is too large';
        flagTerm = 2;
        break;
    end

    if opt.fast
        vk = xk + (cache_yk.Get_ProxGradStep()-xk)/theta;
    end

    xk = cache_yk.Get_ProxGradStep();

    % display stuff

    if opt.display == 1
        Util_PrintProgress(it);
    elseif (opt.display == 2 && mod(it,100) == 0) || opt.display >= 3
        res_curr = norm(cache_yk.Get_FPR(), 'inf')/cache_yk.Get_Gamma();
        obj_curr = cache_yk.Get_FBE();
        fprintf('%6d %7.4e %7.4e %7.4e\n', it, gam, res_curr, obj_curr);
    end

end

time = toc(t0);

if it == opt.maxit
    msgTerm = [msgTerm, 'exceeded maximum iterations'];
    flagTerm = 1;
end

if opt.display == 1
    Util_PrintProgress(it, flagTerm);
elseif opt.display >= 2
    res_curr = norm(cache_yk.Get_FPR(), 'inf')/cache_yk.Get_Gamma();
    obj_curr = cache_yk.Get_FBE();
    fprintf('%6d %7.4e %7.4e %7.4e\n', it, gam, res_curr, obj_curr);
end

% pack up results

out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
out.x = cache_yk.Get_ProxGradStep();
out.iterations = it;
out.operations = ops;
if opt.report
    out.residual = residual(1, 1:it);
    out.objective = objective(1, 1:it);
    out.ts = ts(1, 1:it);
end
out.record = record;
out.gam = gam;
out.adaptive = adaptive;
out.time = time;
