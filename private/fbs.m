function out = fbs(prob, opt)

% initialize output stuff
residual = zeros(1, opt.maxit);
ts = zeros(1, opt.maxit);
objective = zeros(1, opt.maxit);
msgTerm = '';

% display stuff
if opt.display >= 2
    fprintf('%6s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.');
end

if opt.toRecord
    record = [];
end

% initialize operations counter
ops = OpsInit();

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

    cache_yk = CacheInit(prob, yk, gam);
    [cache_yk, ops1] = CacheFBE(cache_yk, gam);
    u = cache_yk.z;
    ops = OpsSum(ops, ops1);

    if it == 1
        cache_0 = cache_yk;
    end

    flagChangedGamma = 0;
    if opt.adaptive || prob.unknownLf
        cache1 = CacheInit(prob, u, gam);
        [cache1, ops1] = CacheEvalf(cache1);
        ops = OpsSum(ops, ops1);
        fz = cache1.fx;
        % increase Lf until a candidate Lipschitz constant is found
        while fz + cache_yk.gz > cache_yk.FBE + 1e-14*(1+abs(cache_yk.FBE))
            flagChangedGamma = 1;
            prob.Lf = prob.Lf*2;
            gam = 1/prob.Lf;
            [cache_yk, ops1] = CacheFBE(cache_yk, gam);
            ops = OpsSum(ops, ops1);
            u = cache_yk.z;
            cache1 = CacheInit(prob, u, gam);
            [cache1, ops1] = CacheEvalf(cache1);
            fz = cache1.fx;
            ops = OpsSum(ops, ops1);
        end
    end

    ts(1, it) = toc(t0);
    residual(1, it) = norm(cache_yk.diff, 'inf')/gam;
    objective(1, it) = cache_yk.FBE;
    if opt.toRecord
        record = [record, opt.record(prob, it, gam, cache_0, cache_yk, ops)];
    end

    % check for termination
    if isnan(cache_yk.normdiff)
        msgTerm = 'something went wrong';
        flagTerm = 1;
        break;
    end
    if residual(1, it) <= 100*eps % cache_yk.normdiff/(1+norm(cache_yk.normdiff)) <= 100*eps
        msgTerm = 'reached optimum (fpr close to eps)';
        flagTerm = 0;
        break;
    end
    if ~flagChangedGamma
        if ~opt.customTerm
            if residual(1, it) <= opt.tol*residual(1, 1)
                msgTerm = 'reached optimum (up to tolerance)';
                flagTerm = 0;
                break;
            end
        else
            flagStop = opt.term(prob, it, gam, cache_0, cache_yk, ops);
            if (prob.unknownLf == 0 || it > 1) && flagStop
                msgTerm = 'reached optimum (custom criterion)';
                flagTerm = 0;
                break;
            end
        end
    end

    xknew = u;

    if opt.fast
        vk = xk + (cache_yk.z-xk)/theta;
    end
    xk = xknew;

    % display stuff
    if opt.display == 1
        if mod(it, 100) == 0
            fprintf('.');
        end
        if mod(it, 4000) == 0
            fprintf('\n');
        end
    elseif opt.display >= 2
        fprintf('%6d %7.4e %7.4e %7.4e\n', it, gam, residual(1,it), objective(1,it));
    end
end

if it == opt.maxit
    flagTerm = 1;
    msgTerm = [msgTerm, 'exceeded maximum iterations'];
end

out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
out.gam = gam;
out.x = cache_yk.z;
out.iterations = it;
out.operations = ops;
out.residual = residual(1, 1:it);
out.objective = objective(1, 1:it);
out.ts = ts(1, 1:it);
out.prob = prob;
out.opt = opt;
if opt.toRecord
    out.record = record;
end
