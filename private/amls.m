function out = amls(prob, opt, varargin)

% initialize operations counter

ops = Ops_Init();

% initialize gamma and sigma

gam = (1-opt.beta)/prob.Lf;

% display header

if opt.display >= 2
    fprintf('%6s%11s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', '||d||', 'tau');
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
        %%% TODO: can we put cache_xbar in place of ~ and recycle it somewhere later?
        [gam, cache_x, ~, ops1, restart] = Cache_BacktrackGamma(cache_x, gam, opt.beta);
        ops = Ops_Sum(ops, ops1);
    else
        [cache_x, ops1] = Cache_ProxGradStep(cache_x, gam);
        ops = Ops_Sum(ops, ops1);
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
    % will be used as threshold for the line-search

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
        msgTerm = 'L is too large';
        flagTerm = 2;
        break;
    end

    % compute search direction and slope

    if it == 1 || restart
        sk = [];
        yk = [];
    end

    [dir_QN, ~, cacheDir] = ...
        opt.methodfun(prob, opt, it, restart, sk, yk, cache_x.FPR, cacheDir);
    dir_FB = -cache_x.FPR;
    
    % FOR DEBUGGING
%     cache_xbar = Cache_Init(prob, cache_x.z, gam);
%     cache_xbar = Cache_FBE(cache_xbar, gam);
%     if cache_xbar.FBE > cache_x.FBE
%         keyboard;
%     end

    % perform line search

    tau = 1.0; % this *must* be 1.0 for this line-search to work
    [cache_x, ops1] = Cache_LineSearch(cache_x, dir_QN);
    ops = Ops_Sum(ops, ops1);
    [cache_w, ops1] = Cache_LineFBE(cache_x, tau, 1);
    ops = Ops_Sum(ops, ops1);
    if opt.adaptive
        [gam, cache_w, cache_wbar, ops1, restart] = Cache_BacktrackGamma(cache_w, gam, opt.beta);
        ops = Ops_Sum(ops, ops1);
    else
        cache_wbar = [];
    end
    if restart, continue; end
    if cache_w.FBE > cache_x.FBE
        [cache_x, ops1] = Cache_LineSearch(cache_x, [], dir_FB);
        ops = Ops_Sum(ops, ops1);
    end
    while cache_w.FBE > cache_x.FBE
        if tau <= 1e-12
            msgTerm = 'line search failed';
            flagTerm = 3;
            break;
        end
        tau = tau/2;
        [cache_w, ops1] = Cache_SegmentFBE(cache_x, tau);
        ops = Ops_Sum(ops, ops1);
        if opt.adaptive
            [gam, cache_w, cache_wbar, ops1, restart] = Cache_BacktrackGamma(cache_w, gam, opt.beta);
            ops = Ops_Sum(ops, ops1);
            if restart, break; end
        end
    end
    if restart, continue; end
    restart = 0;
    if flagTerm
        break;
    end

    % store pair (s, y) to compute next direction

    sk = cache_w.x - cache_x.x;
    yk = cache_w.FPR - cache_x.FPR;

    % update iterate

    if ~isempty(cache_wbar)
        cache_x = cache_wbar;
    else
        cache_x = Cache_Init(prob, cache_w.z, gam);
    end

    % display stuff

    if opt.display == 1
        Util_PrintProgress(it);
    elseif opt.display >= 2
        fprintf('%6d %7.4e %7.4e %7.4e %7.4e %7.4e\n', it, gam, residual(1,it), objective(1,it), norm(dir_QN), tau);
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
if it == opt.maxit
    out.x = cache_x.x;
else
    out.x = cache_x.z;
end
out.iterations = it;
out.operations = ops;
out.residual = residual(1, 1:it);
out.objective = objective(1, 1:it);
out.ts = ts(1, 1:it);
if opt.toRecord, out.record = record; end
out.gam = gam;

end