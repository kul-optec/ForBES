function [cache, ops] = Cache_Evalf(cache)

ops = Ops_Init();

if cache.flagEvalf == 1
    return;
end

prob = cache.prob;

if prob.istheref1
    if prob.isthereC1
        C1x = prob.C1*cache.x;
        ops.C1 = ops.C1 + 1;
        cache.res1x = C1x + prob.d1;
        [cache.f1x, cache.gradf1res1x] = prob.callf1(cache.res1x);
    else
        cache.res1x = cache.x + prob.d1;
        [cache.f1x, cache.gradf1res1x] = prob.callf1(cache.res1x);
        cache.gradf1x = cache.gradf1res1x;
    end
    ops.f1 = ops.f1 + 1;
    ops.gradf1 = ops.gradf1 + 1;
else
    cache.f1x = 0;
end

if prob.istheref2
    if prob.isthereC2
        C2x = prob.C2*cache.x;
        ops.C2 = ops.C2 + 1;
        cache.res2x = C2x + prob.d2;
        f2x = prob.callf2(cache.res2x);
    else
        cache.res2x = cache.x + prob.d2;
        f2x = prob.callf2(cache.res2x);
    end
    ops.f2 = ops.f2 + 1;
    cache.f2x = f2x;
else
    cache.f2x = 0;
end

if prob.istherelin
    cache.flinx = prob.l(:)'*cache.x(:);
    cache.fx = cache.f1x + cache.f2x + cache.flinx;
else
    cache.fx = cache.f1x + cache.f2x;
end

cache.flagEvalf = 1;

end
