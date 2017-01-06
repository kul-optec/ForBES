function out = zerofpr(prob, opt, lsopt)

% initialize operations counter

ops = Ops_Init();

% initialize gamma and sigma

gam = (1-opt.beta)/prob.Lf;
sig = (1-gam*prob.Lf)/(4*gam);

% display header

if opt.display >= 2
    fprintf('%6s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', 'tau');
end

cacheDir.cntSkip = 0;

flagTerm = 0;

MAXIMUM_Lf = 1e14;

t0 = tic();

cache_x = Cache_Init(prob, prob.x0, gam);
restart = 0;

for it = 1:opt.maxit

    % backtracking on gamma

    if opt.adaptive
        [gam, cache_x, cache_xbar, ops1, restart] = Cache_BacktrackGamma(cache_x, gam, opt.beta);
        ops = Ops_Sum(ops, ops1);
    else
        [cache_x, ops1] = Cache_ProxGradStep(cache_x, gam);
        ops = Ops_Sum(ops, ops1);
        cache_xbar = Cache_Init(prob, cache_x.z, gam);
    end

    % trace stuff

    if it == 1
        cache_0 = cache_x;
    end

    ts(1, it) = toc(t0);
    residual(1, it) = norm(cache_x.FPR, 'inf')/cache_x.gam;

    if opt.toRecord
        record(:, it) = opt.record(prob, it, gam, cache_0, cache_x, ops);
    end

    % compute FBE at current point
    % (this should count zero operations)
    % will be used to compute the threshold for the line-search

    [cache_x, ops1] = Cache_FBE(cache_x, gam);
    ops = Ops_Sum(ops, ops1);

    objective(1,it) = cache_x.FBE;

    % check for termination

    if ~restart
        if ~opt.customTerm
            if Cache_StoppingCriterion(cache_x, opt.tol)
                msgTerm = 'reached optimum (up to tolerance)';
                flagTerm = 0;
                break;
            end
        else
            flagStop = opt.term(prob, it, gam, cache_0, cache_x, ops);
            if (opt.adaptive == 0 || it > 1) && flagStop
                msgTerm = 'reached optimum (custom criterion)';
                flagTerm = 0;
                break;
            end
        end
    end
    if prob.Lf >= MAXIMUM_Lf
        msgTerm = ['estimate for Lf became too large: ', num2str(prob.Lf)];
        flagTerm = 1;
        break;
    end

    % compute search direction

    if it == 1 || restart
        sk = [];
        yk = [];
    end

    [cache_xbar, ops1] = Cache_ProxGradStep(cache_xbar, gam);
    ops = Ops_Sum(ops, ops1);

    [dir, tau0, cacheDir] = ...
        opt.methodfun(prob, opt, it, restart, sk, yk, cache_xbar.FPR, cacheDir);
    
    % FOR DEBUGGING
%     cache_xbar = Cache_Init(prob, cache_x.z, gam);
%     cache_xbar = Cache_FBE(cache_xbar, gam);
%     if cache_xbar.FBE > cache_x.FBE
%         keyboard;
%     end

    % perform line search

    ref = cache_x.FBE;
    lin = 0.0;
    const = -sig*cache_x.normFPR^2;
    [tau, cache_tau, ~, ops1, lsopt, ~] = ...
        lsopt.linesearchfun(cache_xbar, dir, 0.0, tau0, lsopt, it, restart, ref, lin, const);
    ops = Ops_Sum(ops, ops1);

    % store pair (s, y) to compute next direction

    sk = cache_tau.x - cache_xbar.x;
    yk = cache_tau.FPR - cache_xbar.FPR;

    % update iterate

    cache_x = cache_tau;

    % display stuff

    if opt.display == 1
        Util_PrintProgress(it);
    elseif opt.display >= 2
        fprintf('%6d %7.4e %7.4e %7.4e %7.4e\n', it, gam, residual(1,it), objective(1,it), tau);
    end

end

if it == opt.maxit
    msgTerm = 'exceeded maximum iterations';
    flagTerm = 1;
end

if opt.display == 1
    Util_PrintProgress(it, flagTerm);
end

% pack up results

out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
out.x = cache_x.z;
out.iterations = it;
out.operations = ops;
out.residual = residual(1, 1:it);
out.objective = objective(1, 1:it);
out.ts = ts(1, 1:it);
if opt.toRecord, out.record = record; end
out.gam = gam;

end
