function out = classical(prob, opt, lsopt)

% initialize operations counter

ops = FBOperations();

% get Lipschitz constant & adaptiveness

[Lf, adaptive] = prob.Get_Lipschitz(opt);

% initialize gamma

gam = (1-opt.beta)/Lf;

% display header

if opt.display >= 2
    fprintf('%6s%11s%11s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', '||dir||', 'slope', 'tau');
end

cacheDir.cntSkip = 0;

t0 = tic();

cache_current = FBCache(prob, prob.x0, gam, ops);
restart = 0;

for it = 1:opt.maxit

    % trace stuff

    if it == 1
        cache_0 = cache_current;
    end

    residual(1, it) = norm(cache_current.Get_FPR(), 'inf')/cache_current.Get_Gamma();
    objective(1, it) = cache_current.Get_FBE();

    ts(1, it) = toc(t0);

    if opt.toRecord
        record(:, it) = opt.record(prob, it, gam, cache_0, cache_current, ops);
    end

    solution = cache_current.z;

    % check for termination

    if isnan(cache_current.normFPR)
        msgTerm = 'something went wrong';
        flagTerm = 1;
        break;
    end
    if ~restart
        if ~opt.customTerm
            if cache_current.Check_StoppingCriterion(opt.tol)
                msgTerm = 'reached optimum (up to tolerance)';
                flagTerm = 0;
                break;
            end
        else
            flagStop = opt.term(prob, it, gam, cache_0, cache_current, ops);
            if (opt.adaptive == 0 || it > 1) && flagStop
                msgTerm = 'reached optimum (custom criterion)';
                flagTerm = 0;
                break;
            end
        end
    end

    % compute pair (s, y) for quasi-Newton updates

    if it > 1
        sk = cache_current.Get_Point() - cache_previous.Get_Point();
        yk = cache_current.Get_GradFBE() - cache_previous.Get_GradFBE();
    else
        sk = [];
        yk = [];
    end

    % compute search direction and slope

    [dir, tau0, cacheDir] = ...
        opt.methodfun(prob, opt, it, restart, sk, yk, cache_current.Get_GradFBE(), cacheDir);
    slope = cache_current.Get_GradFBE()'*dir;

    % perform line search

    [tau, cache_tau, ~, lsopt, flagLS] = ...
        lsopt.linesearchfun(cache_current, dir, slope, tau0, lsopt, adaptive, it, restart);

    % prepare next iteration, store current solution

    restart = 0;
    if flagLS == -1 % gam was too large
        cache_previous = cache_current;
        gam = gam/2;
        restart = 1;
        solution = cache_current.Get_ProxGradStep();
    elseif flagLS > 0 % line-search failed
        flagTerm = 2;
        msgTerm = strcat(['line search failed at iteration ', num2str(it)]);
        break;
    else
        cache_previous = cache_current;
        cache_current = cache_tau;
        solution = cache_tau.z;
    end

    % display stuff

    if opt.display == 1
        Util_PrintProgress(it);
    elseif opt.display >= 2 && mod(it,10) == 0
        fprintf('%6d %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %d\n', it, gam, residual(1,it), objective(1,it), norm(dir), slope, tau, flagLS);
    end

end

if it == opt.maxit
    flagTerm = 1;
    msgTerm = 'exceeded maximum iterations';
end

if opt.display == 1
    Util_PrintProgress(it, flagTerm);
end

% pack up results

out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
out.x = solution;
out.iterations = it;
out.operations = ops;
out.residual = residual(1, 1:it);
out.objective = objective(1, 1:it);
out.ts = ts(1, 1:it);
if opt.toRecord, out.record = record; end
out.gam = gam;
out.skip = cacheDir.cntSkip;
