function out = amls(prob, opt, varargin)

% initialize operations counter

ops = FBOperations();

% get Lipschitz constant & adaptiveness

[Lf, adaptive] = prob.Get_Lipschitz(opt);

% initialize gamma and sigma

gam = (1-opt.beta)/Lf;

% display header

if opt.display >= 2
    fprintf('\n%6s%11s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', '||d||', 'tau');
end

cacheDir.cntSkip = 0;

flagTerm = 0;

MAXIMUM_Lf = 1e14;

t0 = tic();

cache_x = FBCache(prob, prob.x0, gam, ops);
restart = 0;

for it = 1:opt.maxit

    % backtracking on gamma

    if adaptive
        [restart, ~] = cache_x.Backtrack_Gamma(opt.beta);
        gam = cache_x.Get_Gamma();
    end

    % trace stuff

    if it == 1
        cache_0 = cache_x;
    end

    ts(1, it) = toc(t0);
    residual(1, it) = norm(cache_x.Get_FPR(), 'inf')/cache_x.Get_Gamma();
    if opt.toRecord
        record(:, it) = opt.record(prob, it, gam, cache_0, cache_x, ops);
    end

    % compute FBE at current point
    % (this should count zero operations)
    % will be used as threshold for the line-search

    objective(1,it) = cache_x.Get_FBE();

    % check for termination

    if ~restart
        if ~opt.customTerm
            if cache_x.Check_StoppingCriterion(opt.tol)
                msgTerm = 'reached optimum (up to tolerance)';
                flagTerm = 0;
                break;
            end
        else
            flagStop = opt.term(prob, it, gam, cache_0, cache_x, ops);
            if (adaptive == 0 || it > 1) && flagStop
                msgTerm = 'reached optimum (custom criterion)';
                flagTerm = 0;
                break;
            end
        end
    end

    % compute search direction and slope

    if it == 1 || restart
        sk = [];
        yk = [];
    end

    [dir_QN, ~, cacheDir] = ...
        opt.methodfun(prob, opt, it, restart, sk, yk, cache_x.Get_FPR(), cacheDir);
    dir_FB = -cache_x.Get_FPR();

    % perform line search

    tau = 1.0; % this *must* be 1.0 for this line-search to work
    cache_x.Set_Directions(dir_QN);
    cache_w = cache_x.Get_CacheLine(tau, 1);
    if adaptive
        [restart, cache_wbar] = cache_w.Backtrack_Gamma(opt.beta);
        gam = cache_x.Get_Gamma();
    else
        cache_wbar = [];
    end
    if restart, continue; end
    if cache_w.Get_FBE() > cache_x.Get_FBE()
        cache_x.Set_Directions([], dir_FB);
    end
    while cache_w.Get_FBE() > cache_x.Get_FBE()
        if tau <= 1e-14
            msgTerm = 'line search failed';
            flagTerm = 3;
            break;
        end
        tau = tau/2;
        cache_w = cache_x.Get_CacheSegment(tau);
        if adaptive
            [restart, cache_wbar] = cache_w.Backtrack_Gamma(opt.beta);
            gam = cache_x.Get_Gamma();
            if restart, break; end
        end
    end
    if restart, continue; end
    restart = 0;
    if flagTerm
        break;
    end

    % store pair (s, y) to compute next direction

    sk = cache_w.Get_Point() - cache_x.Get_Point();
    yk = cache_w.Get_FPR() - cache_x.Get_FPR();

    % update iterate

    if ~isempty(cache_wbar)
        cache_x = cache_wbar;
    else
        cache_x = FBCache(prob, cache_w.Get_ProxGradStep(), gam, ops);
    end

    % display stuff

    if opt.display == 1
        Util_PrintProgress(it);
    elseif (opt.display == 2 && mod(it,10) == 0) || opt.display >= 3 
        fprintf('%6d %7.4e %7.4e %7.4e %7.4e %7.4e\n', it, gam, residual(1,it), objective(1,it), norm(dir_QN), tau);
    end

end

if it == opt.maxit
    msgTerm = 'exceeded maximum iterations';
    flagTerm = 1;
end

if opt.display == 1
    Util_PrintProgress(it, flagTerm);
elseif opt.display >= 2
    fprintf('%6d %7.4e %7.4e %7.4e\n', it, gam, residual(1,it), objective(1,it));
end

% pack up results

out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
if it == opt.maxit
    out.x = cache_x.Get_Point();
else
    out.x = cache_x.Get_ProxGradStep();
end
out.iterations = it;
out.operations = ops;
out.residual = residual(1, 1:it);
out.objective = objective(1, 1:it);
out.ts = ts(1, 1:it);
if opt.toRecord, out.record = record; end
out.gam = gam;

end
