function out = fbs(prob, opt)
    %% initialize output stuff
    residual = zeros(1, opt.maxit);
    ts = zeros(1, opt.maxit);
    objective = zeros(1, opt.maxit);
    msgTerm = '';
    
    %% display stuff
    if opt.display >= 2
        fprintf('%6s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.');
    end
    
    if opt.toRecord
        record = [];
    end
    
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];
    flagChangedGamma = 0;
    
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
        
        [cache_yk, cnt1] = ForwardBackwardStep(prob, gam, yk);
        u = cache_yk.z;
        cnt = cnt+cnt1;
        
        if it == 1
            cache_0 = cache_yk;
        end
        
        flagChangedGamma = 0;
        if opt.adaptive || prob.unknownLf
            [fz, cnt1] = Evaluatef(prob, u);
            cnt = cnt+cnt1;
            % increase Lf until a candidate Lipschitz constant is found
            while fz + cache_yk.gz > cache_yk.FBE
                flagChangedGamma = 1;
                prob.Lf = prob.Lf*2;
                gam = 1/prob.Lf;
                [cache_yk, cnt1] = ForwardBackwardStep(prob, gam, yk, cache_yk);
                u = cache_yk.z;
                [fz, cnt2] = Evaluatef(prob, cache_yk.z);
                cnt = cnt+cnt1+cnt2;
            end
        end
        
        ts(1, it) = toc(t0);
        residual(1, it) = norm(cache_yk.diff, inf)/gam;
        objective(1, it) = cache_yk.FBE;
        if opt.toRecord
            record = [record, opt.record(prob, it, gam, cache_0, cache_yk, cnt)];
        end
        
        %% check for termination
        if ~flagChangedGamma
            if ~opt.customTerm
                % From sec. 8.2.3.2 of Gill, Murray, Wright (1982).
                absFBE = abs(cache_current.FBE);
                flagStop = residual(1, it) <= nthroot(opt.tol, 3)*(1+absFBE) && ...
                    norm(cache_yk1.z-cache_yk.z, inf) < sqrt(opt.tol)*(1+norm(cache_yk.z, inf)) && ...
                    abs(objective(1, it-1)-objective(1, it)) < opt.tol*(1+absFBE);
                if residual(1, it) <= 10*sqrt(eps) || ((prob.unknownLf == 0 || it > 1) && flagStop)
                    msgTerm = 'reached optimum (up to tolerance)';
                    flagTerm = 0;
                    break;
                end
            else
                flagStop = opt.term(prob, it, gam, cache_0, cache_yk, cnt);
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
        cache_yk1 = cache_yk;
        
        %% display stuff
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
    out.operations.cnt_f1 = cnt(1);
    out.operations.cnt_C1 = cnt(2);
    out.operations.cnt_f2 = cnt(4);
    out.operations.cnt_C2 = cnt(3);
    out.operations.cnt_g = cnt(5);
    out.residual = residual(1, 1:it);
    out.objective = objective(1, 1:it);
    out.ts = ts(1, 1:it);
    out.prob = prob;
    out.opt = opt;
    if opt.toRecord
        out.record = record;
    end
end

%% compute the FBE value and store reusable quantities
function [cache, cnt] = ForwardBackwardStep(prob, gam, x, cache)
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];
    f1x = 0; gradf1x = 0;
    f2x = 0; gradf2x = 0;
    cache.gam = gam;
    if nargin < 4
        cache.x = x;
        if prob.istheref1
            if prob.isthereC1
                if prob.isC1fun, C1x = prob.C1(cache.x);
                else C1x = prob.C1*cache.x; end
                cache.res1x = C1x - prob.d1;
                if prob.isQfun, cache.Qres1x = prob.Q(cache.res1x);
                else cache.Qres1x = prob.Q*cache.res1x; end
                if prob.isC1fun, gradf1x = prob.C1t(cache.Qres1x + prob.q);
                else gradf1x = prob.C1'*(cache.Qres1x + prob.q); end
                cnt(2) = cnt(2)+2;
            else
                cache.res1x = cache.x - prob.d1;
                if prob.isQfun, cache.Qres1x = prob.Q(cache.res1x);
                else cache.Qres1x = prob.Q*cache.res1x; end
                gradf1x = cache.Qres1x + prob.q;
            end
            cnt(1) = cnt(1)+1;
            cache.gradf1x = gradf1x;
            f1x = 0.5*(cache.res1x'*cache.Qres1x) + prob.q'*cache.res1x;
            cache.f1x = f1x;
        end
        if prob.istheref2
            if prob.isthereC2
                if prob.isC2fun, C2x = prob.C2(cache.x);
                else C2x = prob.C2*cache.x; end
                cache.res2x = C2x - prob.d2;
%                 if prob.useHessian
%                     [f2x, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
%                 else
                    [f2x, gradf2res2x] = prob.callf2(cache.res2x);
%                 end
                if prob.isC2fun, gradf2x = prob.C2t(gradf2res2x);
                else gradf2x = prob.C2'*gradf2res2x; end
                cnt(3) = cnt(3)+2;
            else
                cache.res2x = cache.x - prob.d2;
%                 if prob.useHessian
%                     [f2x, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
%                 else
                    [f2x, gradf2res2x] = prob.callf2(cache.res2x);
%                 end
                gradf2x = gradf2res2x;
            end
            cnt(4) = cnt(4)+1;
            cache.gradf2x = gradf2x;
        end
        if prob.istherelin
            cache.flinx = prob.l'*x;
            cache.fx = f1x + f2x + cache.flinx;
            cache.gradfx = gradf1x + gradf2x + prob.l;
        else
            cache.fx = f1x + f2x;
            cache.gradfx = gradf1x + gradf2x;
        end
    end
    y = cache.x - gam*cache.gradfx;
    [cache.z, cache.gz] = prob.callg(y, gam);
    cnt(5) = cnt(5)+1;
    cache.diff = cache.z-cache.x;
    sqnormdiff = cache.diff'*cache.diff;
    cache.normdiff = sqrt(sqnormdiff);
    cache.FBE = cache.fx + cache.gz + ...
        cache.gradfx'*cache.diff + ...
        (0.5/gam)*sqnormdiff;
end

