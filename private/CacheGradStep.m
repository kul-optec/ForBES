function [cache, ops] = CacheGradStep(cache, gam)

ops = OpsInit();

if nargin < 2
    gam = cache.gam;
end

if cache.flagGradStep == 1
    if cache.gam ~= gam
        cache.gam = gam;
        cache.y = cache.x - gam*cache.gradfx;
        cache.flagProxGradStep = 0;
        cache.flagFBE = 0;
        cache.flagGradFBE = 0;
    end
    return;
end

if cache.flagEvalf == 0
    [cache, ops] = CacheEvalf(cache);
end

cache.gam = gam;
prob = cache.prob;

if prob.istheref1
    if prob.isthereC1
        if prob.isC1fun, cache.gradf1x = prob.C1t(cache.gradf1res1x);
        else cache.gradf1x = prob.C1'*(cache.gradf1res1x); end
        ops.C1 = ops.C1 + 1;
    end
else
    cache.gradf1x = 0;
end

if prob.istheref2
    if prob.isthereC2
        if prob.useHessian
            [~, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
        else
            [~, gradf2res2x] = prob.callf2(cache.res2x);
            cache.gradf2res2x = gradf2res2x;
        end
        if prob.isC2fun, cache.gradf2x = prob.C2t(gradf2res2x);
        else cache.gradf2x = prob.C2'*gradf2res2x; end
        ops.C2 = ops.C2 + 1;
    else
        if prob.useHessian
            [~, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
        else
            [~, gradf2res2x] = prob.callf2(cache.res2x);
            cache.gradf2res2x = gradf2res2x;
        end
        cache.gradf2x = gradf2res2x;
    end
    ops.gradf2 = ops.gradf2 + 1;
else
    cache.gradf2x = 0;
end

if prob.istherelin
    cache.gradfx = cache.gradf1x + cache.gradf2x + prob.l;
else
    cache.gradfx = cache.gradf1x + cache.gradf2x;
end

cache.y = cache.x - gam*cache.gradfx;

cache.flagGradStep = 1;
