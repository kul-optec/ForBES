function [cache, ops] = CacheGradStep(cache, gam)

ops = OpsInit();

if cache.flagEvalf == 0
    [cache, ops] = CacheEvalf(cache);
end

% if the gradient step was performed already, then compute it only
% if gam (the stepsize) has changed, otherwise don't do anything
if cache.flagGradStep == 1
    if cache.gam ~= gam
        cache.gam = gam;
        cache.y = cache.x - gam*cache.gradfx;
    end
    return;
end

cache.gam = gam;
prob = cache.prob;
gradf1x = 0;
gradf2x = 0;

if prob.istheref1
    if prob.isthereC1
        if prob.isC1fun, gradf1x = prob.C1t(cache.Qres1x + prob.q);
        else gradf1x = prob.C1'*(cache.Qres1x + prob.q); end
        ops.C1 = ops.C1 + 1;
    else
        gradf1x = cache.Qres1x + prob.q;
    end
    cache.gradf1x = gradf1x;
end

if prob.istheref2
    if prob.isthereC2
        if prob.useHessian
            [~, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
        else
            [~, gradf2res2x] = prob.callf2(cache.res2x);
            cache.gradf2res2x = gradf2res2x;
        end
        if prob.isC2fun, gradf2x = prob.C2t(gradf2res2x);
        else gradf2x = prob.C2'*gradf2res2x; end
        ops.C2 = ops.C2 + 1;
    else
        if prob.useHessian
            [~, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
        else
            [~, gradf2res2x] = prob.callf2(cache.res2x);
            cache.gradf2res2x = gradf2res2x;
        end
        gradf2x = gradf2res2x;
    end
    ops.gradf2 = ops.gradf2 + 1;
    cache.gradf2x = gradf2x;
end

if prob.istherelin
    cache.gradfx = gradf1x + gradf2x + prob.l;
else
    cache.gradfx = gradf1x + gradf2x;
end

cache.y = cache.x - gam*cache.gradfx;

cache.flagGradStep = 1;
