function out = classical(prob, opt, lsopt)

% initialize operations counter

ops = Ops_Init();

% initialize gamma

gam = (1-opt.beta)/prob.Lf;

% display header

if opt.display >= 2
    fprintf('%6s%11s%11s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', '||dir||', 'slope', 'tau');
end

cacheDir.cntSkip = 0;
hasGammaChanged = 0;

cache_current = Cache_Init(prob, prob.x0, gam);

t0 = tic();

for it = 1:opt.maxit

    % compute FBE

    [cache_current, ops1] = Cache_FBE(cache_current, gam);
    ops = Ops_Sum(ops, ops1);

    % store initial cache

    if it == 1
        cache_0 = cache_current;
    end

    % trace stuff

    ts(1, it) = toc(t0);
    residual(1, it) = norm(cache_current.FPR, 'inf')/gam;
    objective(1, it) = cache_current.FBE;
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
    if ~hasGammaChanged
        if ~opt.customTerm
            if Cache_StoppingCriterion(cache_current, opt.tol)
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

    % compute gradient of the FBE

    [cache_current, ops1] = Cache_GradFBE(cache_current, gam);
    ops = Ops_Sum(ops, ops1);

    % compute pair (s, y) for quasi-Newton updates

    if it > 1
        sk = cache_current.x - cache_previous.x;
        yk = cache_current.gradFBE - cache_previous.gradFBE;
    else
        sk = [];
        yk = [];
    end

    % compute search direction and slope

    [dir, tau0, cacheDir] = ...
        opt.methodfun(prob, opt, it, hasGammaChanged, sk, yk, cache_current.gradFBE, cacheDir);
    slope = cache_current.gradFBE'*dir;

    % perform line search

    [tau, cache_tau, ~, ops1, lsopt, flagLS] = ...
        lsopt.linesearchfun(cache_current, dir, slope, tau0, lsopt, it, hasGammaChanged);
    ops = Ops_Sum(ops, ops1);

    % prepare next iteration, store current solution

    hasGammaChanged = 0;
    if flagLS == -1 % gam was too large
        cache_previous = cache_current;
        prob.Lf = prob.Lf*2; gam = gam/2;
        hasGammaChanged = 1;
        solution = cache_current.z;
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
    elseif opt.display >= 2
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
