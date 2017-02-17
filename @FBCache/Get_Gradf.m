% GET_GRADF(cache)
%
%   Returns the gradient of f, the smooth term of the problem, at the point x
%   to which cache refers.

function gradfx = Get_Gradf(cache)

if cache.flagGradf == true
    gradfx = cache.gradfx;
    return;
end

if cache.flagEvalf == false
    cache.Get_f();
end

prob = cache.prob;

if prob.istheref1
    if prob.isthereC1
        cache.gradf1x = prob.C1'*cache.gradf1res1x;
        if cache.flagOps, cache.ops.addC1(); end
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
        if cache.flagOps, cache.ops.addC2(); end
    else
        if prob.useHessian
            [~, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
        else
            [~, gradf2res2x] = prob.callf2(cache.res2x);
            cache.gradf2res2x = gradf2res2x;
        end
        cache.gradf2x = gradf2res2x;
    end
    if cache.flagOps, cache.ops.addgradf2(); end
else
    cache.gradf2x = 0.0;
end

if prob.istherelin
    cache.gradfx = cache.gradf1x + cache.gradf2x + prob.lin;
else
    cache.gradfx = cache.gradf1x + cache.gradf2x;
end

cache.flagGradf = true;
gradfx = cache.gradfx;
