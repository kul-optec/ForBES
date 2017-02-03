function out = minfbe(prob, opt, lsopt)

% initialize operations counter

ops = FBOperations();

% get Lipschitz constant & adaptiveness

[Lf, adaptive] = prob.Get_Lipschitz(opt);

% initialize gamma

gam = (1-opt.beta)/Lf;

% display header

if opt.display >= 2
    fprintf('\n%6s%11s%11s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', '||dir||', 'slope', 'tau');
end

cache_dir.cntSkip = 0;
restart = 0;

cache_current = FBCache(prob, prob.x0, gam, ops);

t0 = tic();

for it = 1:opt.maxit

    % store initial cache

    if it == 1
        cache_0 = cache_current;
    end

    % trace stuff

    objective(1, it) = cache_current.Get_FBE(); %cache_current.FBE;
    ts(1, it) = toc(t0);
    residual(1, it) = norm(cache_current.Get_FPR(), 'inf')/cache_current.Get_Gamma();
    if opt.toRecord
        record(:, it) = opt.record(prob, it, gam, cache_0, cache_current, ops);
    end

    solution = cache_current.Get_ProxGradStep();

    % check for termination

    if ~restart
        if ~opt.customTerm
            if cache_current.Check_StoppingCriterion(opt.tol)
                msgTerm = 'reached optimum (up to tolerance)';
                flagTerm = 0;
                break;
            end
        else
            flagStop = opt.term(prob, it, gam, cache_0, cache_current, ops);
            if (adaptive == 0 || it > 1) && flagStop
                msgTerm = 'reached optimum (custom criterion)';
                flagTerm = 0;
                break;
            end
        end
    end

    % store pair (s, y) to compute direction

    if it > 1
        if opt.memopt == 1
            sk = cache_current.Get_Point() - cache_previous.Get_Point();
            yk = cache_current.Get_GradFBE() - cache_previous.Get_GradFBE();
        elseif opt.memopt == 2
            sk = cache_tau.Get_Point() - cache_previous.Get_Point();
            yk = cache_tau.Get_GradFBE() - cache_previous.Get_GradFBE();
        end
    else
        sk = [];
        yk = [];
    end

    % compute search direction and slope

    [dir, tau0, cache_dir] = ...
        opt.methodfun(prob, opt, it, restart, sk, yk, cache_current.Get_GradFBE(), cache_dir);
    slope = cache_current.Get_GradFBE()'*dir;

    % perform line search

    [tau, cache_tau, cache_tau1, lsopt, flagLS] = ...
        lsopt.linesearchfun(cache_current, dir, slope, tau0, lsopt, adaptive, it, restart);

    % prepare next iteration, store current solution

    restart = 0;
    if flagLS == -1 % gam was too large
        cache_previous = cache_current;
        gam = gam/2;
        restart = 1;
        solution = cache_current.Get_ProxGradStep();
    elseif flagLS > 0 % line-search failed
        cache_previous = cache_current;
        cache_current = FBCache(prob, cache_current.Get_ProxGradStep(), gam, ops);
        solution = cache_current.Get_Point();
    else
        cache_previous = cache_current;
        if ~isempty(cache_tau1)
            solution = cache_current.Get_ProxGradStep();
            cache_current = cache_tau1;
        else
            solution = cache_tau.Get_ProxGradStep();
            cache_current = FBCache(prob, cache_tau.Get_ProxGradStep(), gam, ops);
        end
    end

    % display stuff

    if opt.display == 1
        Util_PrintProgress(it);
    elseif (opt.display == 2 && mod(it,10) == 0) || opt.display >= 3 
        fprintf('%6d %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %d\n', it, gam, residual(1,it), objective(1,it), norm(dir), slope, tau, flagLS);
    end

end

if it == opt.maxit
    flagTerm = 1;
    msgTerm = 'exceeded maximum iterations';
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
out.x = solution;
out.iterations = it;
out.operations = ops;
out.residual = residual(1, 1:it);
out.objective = objective(1, 1:it);
out.ts = ts(1, 1:it);
if opt.toRecord, out.record = record; end
out.gam = gam;
out.skip = cache_dir.cntSkip;
