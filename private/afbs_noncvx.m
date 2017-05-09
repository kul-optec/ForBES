function out = afbs_noncvx(prob, opt, varargin)

MAXIMUM_Lf = 1e14;

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

gam_x = (1-opt.beta)/Lf;
gam_y = (1-opt.beta)/Lf;

t1 = 1; t0 = 0;
eta = 0.5; % eta = 0.0 for monotone, > 0.0 for nonmonotone
q = 1.0;
del = 1e-2; % the smaller this is, the more extrapolations are accepted
c = opt.F0;

cache_z1 = FBCache(prob, prob.x0, gam_x, ops);
cache_x1 = FBCache(prob, prob.x0, gam_x, ops);
cache_x0 = FBCache(prob, prob.x0, gam_x, ops);

tic0 = tic();

for it = 1:opt.maxit
    
    z1 = cache_z1.Get_Point();
    x1 = cache_x1.Get_Point();
    x0 = cache_x0.Get_Point();
    
    y1 = x1 + t0/t1*(z1 - x1) + (t0-1)/t1*(x1 - x0);
    cache_y1 = FBCache(prob, y1, gam_y, ops);
    
    if adaptive
        [hasGammaChanged_y, ~] = cache_y1.Backtrack_Gamma(opt.beta);
        gam_y = cache_y1.Get_Gamma();
    end
    
    z1 = cache_y1.Get_ProxGradStep();
    cache_z1 = FBCache(prob, z1, gam_x, ops);
    
    % record values, stopping criterion
    
    if it == 1
        cache_0 = cache_y1;
    end
    
    if opt.toRecord
        record = [record, opt.record(prob, it, cache_0, cache_y1)];
    end
    
    if opt.report
        objective(1, it) = cache_y1.Get_FBE();
        residual(1, it) = norm(cache_y1.Get_FPR(), 'inf')/cache_y1.Get_Gamma();
        ts(1, it) = toc(tic0);
    end
    
    if cache_y1.Check_StoppingCriterion(opt.tol)
        msgTerm = 'reached optimum (up to tolerance)';
        flagTerm = 0;
        break;
    end
    
    %%%
    
    F_z1 = cache_y1.Get_g() + cache_z1.Get_f();
    
    cache_x0 = cache_x1;
    
    if F_z1 <= c - del*cache_y1.Get_NormFPR()^2
        % extrapolation is accepted
        cache_x1 = cache_z1;
        F_x1 = F_z1;
    else
        if adaptive
            [hasGammaChanged_x, ~] = cache_x1.Backtrack_Gamma(opt.beta);
            gam_x = cache_x1.Get_Gamma();
            if hasGammaChanged_x
                cache_z1.Set_Gamma(gam_x);
            end
        end
        v1 = cache_x1.Get_ProxGradStep();
        cache_v1 = FBCache(prob, v1, gam_x, ops);
        F_v1 = cache_v1.Get_f() + cache_x1.Get_g();
        if F_z1 <= F_v1
            % extrapolation is accepted
            cache_x1 = cache_z1;
            F_x1 = F_z1;
        else
            % extrapolation is rejected (ordinary FB step)
            cache_x1 = cache_v1;
            F_x1 = F_v1;
        end
    end
    
    t0 = t1;
    t1 = (sqrt(4*(t1^2) + 1) + 1)/2;
    q1 = eta*q + 1;
    c = (eta*q*c + F_x1)/q1;
    q = q1;
    
    % display stuff

    if opt.display == 1
        Util_PrintProgress(it);
    elseif (opt.display == 2 && mod(it,100) == 0) || opt.display >= 3
        res_curr = norm(cache_y1.Get_FPR(), 'inf')/cache_y1.Get_Gamma();
        obj_curr = cache_y1.Get_FBE();
        fprintf('%6d %7.4e %7.4e %7.4e\n', it, cache_y1.Get_Gamma(), res_curr, obj_curr);
    end

end

time = toc(tic0);

if it == opt.maxit
    msgTerm = [msgTerm, 'exceeded maximum iterations'];
    flagTerm = 1;
end

if opt.display == 1
    Util_PrintProgress(it, flagTerm);
elseif opt.display >= 2
    res_curr = norm(cache_y1.Get_FPR(), 'inf')/cache_y1.Get_Gamma();
    obj_curr = cache_y1.Get_FBE();
    fprintf('%6d %7.4e %7.4e %7.4e\n', it, cache_y1.Get_Gamma(), res_curr, obj_curr);
end

% pack up results

out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
out.x = cache_y1.Get_ProxGradStep();
out.iterations = it;
out.operations = ops;
if opt.report
    out.residual = residual(1, 1:it);
    out.objective = objective(1, 1:it);
    out.ts = ts(1, 1:it);
end
out.record = record;
out.gam = cache_y1.Get_Gamma();
out.adaptive = adaptive;
out.time = time;
