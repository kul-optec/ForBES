function out = zerofprChainBFGS(prob, opt, lsopt)

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

% get Lipschitz constant & adaptiveness

[Lf, adaptive] = prob.Get_Lipschitz(opt);

% initialize gamma and sigma

gam = (1-opt.beta)/Lf;
sig = opt.beta/(4*gam);

% display header

if opt.display >= 2
    fprintf('\n%s', opt.name);
    fprintf('\n%6s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', 'tau');
end

cacheDir.cntSkip = 0;

flagTerm = 0;

t0 = tic();

cache_x = FBCache(prob, prob.x0, gam, ops);
restart = 0;

for it = 1:opt.maxit

    % backtracking on gamma

    if adaptive
        [restart, cache_xbar] = cache_x.Backtrack_Gamma(opt.beta);
        gam = cache_x.Get_Gamma();
        sig = opt.beta/(4*gam);
    else
        x_bar = cache_x.Get_ProxGradStep();
        cache_xbar = FBCache(prob, x_bar, cache_x.Get_Gamma(), ops);
    end

    if opt.report
        objective(1,it) = cache_x.Get_FBE();
        residual(1, it) = norm(cache_x.Get_FPR(), 'inf')/cache_x.Get_Gamma();
        ts(1, it) = toc(t0);
    end

    if it == 1
        cache_0 = cache_x;
    end

    if opt.toRecord
        record(:, it) = opt.record(prob, it, cache_0, cache_x);
    end

    % check for termination

    if ~restart
        if ~opt.customTerm
            if cache_x.Check_StoppingCriterion(opt.tol)
                msgTerm = 'reached optimum (up to tolerance)';
                flagTerm = 0;
                break;
            end
        else
            flagStop = opt.term(prob, it, cache_0, cache_x);
            if (adaptive == 0 || it > 1) && flagStop
                msgTerm = 'reached optimum (custom criterion)';
                flagTerm = 0;
                break;
            end
        end
    end

    % compute search direction

    if it == 1 || restart
      B = 1;
    else
      % next check is to initialize the matrices B1 and B2
      % this initialization should be performed also when restart = 1
      % but fixed gam it is OK like this
      if it == 2
        B1 = (y1'*s1)/(y1'*y1);
        B2 = (y2'*s2)/(y2'*y2);
      end

      % monotonicity check
      YS1 = y1'*s1;
      YS2 = y2'*s2;

      if YS1 > 0 && YS2 > 0
        B1s1 = B1*s1;
        s1B1s1 = s1'*B1s1;
        B1 = B1 + (y1*y1')/YS1 - (B1s1*B1s1')/s1B1s1;

        B2s2 = B2*s2;
        s2B2s2 = s2'*B2s2;
        B2 = B2 + (y2*y2')/YS2 - (B2s2*B2s2')/s2B2s2;
      else
        cacheDir.cntSkip = cacheDir.cntSkip + 1;
      end

      B = eye(prod(prob.n))-B2*B1;
    end

    direction = -B\cache_xbar.Get_FPR();
    tau0 = 1;

    % perform line search

    ref = cache_x.Get_FBE();
    lin = 0.0;
    const = -sig*cache_x.Get_NormFPR()^2;
    [tau, cache_tau, ~, lsopt, ~] = ...
        lsopt.linesearchfun(cache_xbar, direction, 0.0, tau0, lsopt, false, it, restart, ref, lin, const);

    % store pair (s, y) to compute next direction

    s1 = cache_tau.Get_Point() - cache_xbar.Get_Point();
    y1 = cache_tau.Get_GradStep() - cache_xbar.Get_GradStep();
    s2 = y1;
    y2 = cache_tau.Get_ProxGradStep() - cache_xbar.Get_ProxGradStep();

    % update iterate

    cache_x = cache_tau;

    % display stuff

    if opt.display == 1
        Util_PrintProgress(it);
    elseif (opt.display == 2 && mod(it,10) == 0) || opt.display >= 3
        res_curr = norm(cache_x.Get_FPR(), 'inf')/cache_x.Get_Gamma();
        obj_curr = cache_x.Get_FBE();
        fprintf('%6d %7.4e %7.4e %7.4e %7.4e\n', it, gam, res_curr, obj_curr, tau);
    end

end

time = toc(t0);

if it == opt.maxit
    msgTerm = 'exceeded maximum iterations';
    flagTerm = 1;
end

if opt.display == 1
    Util_PrintProgress(it, flagTerm);
elseif opt.display >= 2
    res_curr = norm(cache_x.Get_FPR(), 'inf')/cache_x.Get_Gamma();
    obj_curr = cache_x.Get_FBE();
    fprintf('%6d %7.4e %7.4e %7.4e\n', it, gam, res_curr, obj_curr);
end

% pack up results

out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
out.x = cache_x.Get_ProxGradStep();
out.iterations = it;
out.operations = ops;
if opt.report
    out.residual = residual(1, 1:it);
    out.objective = objective(1, 1:it);
    out.ts = ts(1, 1:it);
end
if opt.toRecord, out.record = record; end
out.gam = gam;
out.adaptive = adaptive;
out.time = time;
out.cacheDir = cacheDir;
