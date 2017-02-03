function y = Get_GradStep(cache, gam)

if nargin < 2
    gam = cache.gam;
end

if cache.flagGradStep
    if cache.gam ~= gam
        cache.gam = gam;
        cache.y = cache.x - gam*cache.gradfx;
        cache.flagProxGradStep = 0;
        cache.flagFBE = 0;
        cache.flagGradFBE = 0;
    end
    y = cache.y;
    return;
end

if cache.flagEvalf == 0
    cache.Get_f();
end

cache.gam = gam;
prob = cache.prob;

if prob.istheref1
    if prob.isthereC1
        cache.gradf1x = prob.C1'*cache.gradf1res1x;
        cache.ops.addC1();
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
        cache.ops.addC2();
    else
        if prob.useHessian
            [~, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
        else
            [~, gradf2res2x] = prob.callf2(cache.res2x);
            cache.gradf2res2x = gradf2res2x;
        end
        cache.gradf2x = gradf2res2x;
    end
    cache.ops.addgradf2();
else
    cache.gradf2x = 0.0;
end

if prob.istherelin
    cache.gradfx = cache.gradf1x + cache.gradf2x + prob.lin;
else
    cache.gradfx = cache.gradf1x + cache.gradf2x;
end

cache.y = cache.x - gam*cache.gradfx;
cache.flagGradStep = true;
y = cache.y;
