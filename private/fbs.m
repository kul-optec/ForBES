function out = fbs(prob, opt, varargin)

MAXIMUM_Lf = 1e14;

% initialize output stuff

residual = zeros(1, opt.maxit);
ts = zeros(1, opt.maxit);
objective = zeros(1, opt.maxit);
msgTerm = '';
record = [];

% display stuff

if opt.display >= 2
    fprintf('%6s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.');
end

% initialize operations counter

ops = Ops_Init();

% gam = (1-opt.beta)/prob.Lf;
gam = 1/prob.Lf;
xk = prob.x0;
vk = prob.x0;

t0 = tic();

for it = 1:opt.maxit

    if opt.fast
        if prob.muf == 0
            theta = 2/(it+1); % since it starts from 1
        else
            theta = sqrt(prob.muf/prob.Lf);
        end
        yk = (1-theta)*xk+theta*vk;
    else
        yk = xk;
    end

    cache_yk = Cache_Init(prob, yk, gam);

    if it == 1
        cache_0 = cache_yk;
    end

    hasGammaChanged = 0;
    if opt.adaptive
        [isGammaOK, cache_yk, ~, ops1] = Cache_CheckGamma(cache_yk, gam, opt.beta);
        ops = Ops_Sum(ops, ops1);
        while ~isGammaOK
            prob.Lf = 2*prob.Lf; gam = gam/2;
            hasGammaChanged = 1;
            [isGammaOK, cache_yk, ~, ops1] = Cache_CheckGamma(cache_yk, gam, opt.beta);
            ops = Ops_Sum(ops, ops1);
        end
    end

    [cache_yk, ops1] = Cache_FBE(cache_yk, gam);
    ops = Ops_Sum(ops, ops1);

    ts(1, it) = toc(t0);
    residual(1, it) = norm(cache_yk.FPR, 'inf')/cache_yk.gam;
    objective(1, it) = cache_yk.FBE;
    if opt.toRecord
        record = [record, opt.record(prob, it, gam, cache_0, cache_yk, ops)];
    end

    if ~hasGammaChanged
        if ~opt.customTerm
            if Cache_StoppingCriterion(cache_yk, opt.tol)
                msgTerm = 'reached optimum (up to tolerance)';
                flagTerm = 0;
                break;
            end
        else
            flagStop = opt.term(prob, it, gam, cache_0, cache_yk, ops);
            if (opt.adaptive == 0 || it > 1) && flagStop
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
        vk = xk + (cache_yk.z-xk)/theta;
    end

    xk = cache_yk.z;

    % display stuff

    if opt.display == 1
        Util_PrintProgress(it);
    elseif opt.display >= 2
        fprintf('%6d %7.4e %7.4e %7.4e\n', it, gam, residual(1,it), objective(1,it));
    end

end

if it == opt.maxit
    msgTerm = [msgTerm, 'exceeded maximum iterations'];
    flagTerm = 1;
end

if opt.display == 1
    Util_PrintProgress(it, flagTerm);
end

% pack up results

out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
out.x = cache_yk.z;
out.iterations = it;
out.operations = ops;
out.residual = residual(1, 1:it);
out.objective = objective(1, 1:it);
out.ts = ts(1, 1:it);
out.record = record;
out.gam = gam;

end
