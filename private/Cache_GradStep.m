function [cache, ops] = Cache_GradStep(cache, gam)

ops = Ops_Init();

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
    [cache, ops] = Cache_Evalf(cache);
end

cache.gam = gam;
prob = cache.prob;

if prob.istheref1
    if prob.isthereC1
        cache.gradf1x = prob.C1'*cache.gradf1res1x;
        ops.C1 = ops.C1 + 1;
    end
else
    cache.gradf1x = 0.0;
end

if prob.istheref2
    if prob.isthereC2
        if prob.useHessian
            [~, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
        else
            [~, gradf2res2x] = prob.callf2(cache.res2x);
            cache.gradf2res2x = gradf2res2x;
        end
        cache.gradf2x = prob.C2'*gradf2res2x;
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
    cache.gradf2x = 0.0;
end

if prob.istherelin
    cache.gradfx = cache.gradf1x + cache.gradf2x + prob.l;
else
    cache.gradfx = cache.gradf1x + cache.gradf2x;
end

cache.y = cache.x - gam*cache.gradfx;

cache.flagGradStep = 1;

end
