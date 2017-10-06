% GET_F(cache)
%
%   Returns the value of f, the smooth term of the problem, evaluated at x, the
%   point to which cache refers.

function fx = Get_f(cache)

if cache.flagEvalf
    fx = cache.fx;
    return;
end

prob = cache.prob;

if prob.istheref1
    if prob.isthereC1
        C1x = prob.C1*cache.x;
        if cache.flagOps, cache.ops.addC1(); end
        cache.res1x = C1x + prob.d1;
        [cache.gradf1res1x, cache.f1x] = prob.f1.gradient(cache.res1x);
    else
        cache.res1x = cache.x + prob.d1;
        [cache.gradf1res1x, cache.f1x] = prob.f1.gradient(cache.res1x);
        cache.gradf1x = cache.gradf1res1x;
    end
    if cache.flagOps
        cache.ops.addf1();
        cache.ops.addgradf1();
    end
else
    cache.f1x = 0;
end

if prob.istheref2
    if prob.isthereC2
        C2x = prob.C2*cache.x;
        if cache.flagOps, cache.ops.addC2(); end
        cache.res2x = C2x + prob.d2;
        [~, f2x] = prob.f2.gradient(cache.res2x);
    else
        cache.res2x = cache.x + prob.d2;
        [~, f2x] = prob.f2.gradient(cache.res2x);
    end
    if cache.flagOps, cache.ops.addf2(); end
    cache.f2x = f2x;
else
    cache.f2x = 0;
end

if prob.istherelin
    cache.flinx = prob.lin(:)'*cache.x(:);
    cache.fx = cache.f1x + cache.f2x + cache.flinx;
else
    cache.fx = cache.f1x + cache.f2x;
end

cache.flagEvalf = true;
fx = cache.fx;
