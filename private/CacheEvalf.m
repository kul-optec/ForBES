function [cache, ops] = CacheEvalf(cache)

ops = OpsInit();

if cache.flagEvalf == 1
    return;
end

prob = cache.prob;
f1x = 0;
f2x = 0;

if prob.istheref1
    if prob.isthereC1
        if prob.isC1fun, C1x = prob.C1(cache.x);
        else C1x = prob.C1*cache.x; end
        cache.res1x = C1x - prob.d1;
        if prob.isQfun, cache.Qres1x = prob.Q(cache.res1x);
        else cache.Qres1x = prob.Q*cache.res1x; end
        ops.C1 = ops.C1 + 1;
    else
        cache.res1x = cache.x - prob.d1;
        if prob.isQfun, cache.Qres1x = prob.Q(cache.res1x);
        else cache.Qres1x = prob.Q*cache.res1x; end
    end
    ops.Q = ops.Q + 1;
    f1x = 0.5*(cache.res1x'*cache.Qres1x) + prob.q'*cache.res1x;
    cache.f1x = f1x;
end

if prob.istheref2
    if prob.isthereC2
        if prob.isC2fun, C2x = prob.C2(cache.x);
        else C2x = prob.C2*cache.x; end
        cache.res2x = C2x - prob.d2;
        f2x = prob.callf2(cache.res2x);
        ops.C2 = ops.C2 + 1;
    else
        cache.res2x = cache.x - prob.d2;
        f2x = prob.callf2(cache.res2x);
    end
    ops.f2 = ops.f2 + 1;
    cache.f2x = f2x;
end

if prob.istherelin
    cache.flinx = prob.l'*x;
    cache.fx = f1x + f2x + cache.flinx;
else
    cache.fx = f1x + f2x;
end

cache.flagEvalf = 1;
