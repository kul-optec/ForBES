function [cache, ops] = CacheLineSearch(cache, dir)

if cache.flagLineSearch == 1
    if norm(dir-cache.dir) == 0
        ops = OpsInit();
        return;
    end
end

if cache.flagGradStep == 0
    [cache, ops] = CacheGradStep(cache, cache.gam);
else
    ops = OpsInit();
end

prob = cache.prob;
cache.dir = dir;

if prob.istheref1
    if prob.isthereC1
        if prob.isC1fun, cache.C1dir = prob.C1(dir);
        else cache.C1dir = prob.C1*dir; end
        if prob.isQfun, cache.QC1dir = prob.Q(cache.C1dir);
        else cache.QC1dir = prob.Q*cache.C1dir; end
        if prob.isC1fun, cache.C1tQC1dir = prob.C1t(cache.QC1dir);
        else cache.C1tQC1dir = prob.C1'*cache.QC1dir; end
        ops.C1 = ops.C1 + 2;
    else
        cache.C1dir = dir;
        if prob.isQfun, cache.QC1dir = prob.Q(cache.C1dir);
        else cache.QC1dir = prob.Q*cache.C1dir; end
        cache.C1tQC1dir = cache.QC1dir;
    end
    ops.Q = ops.Q + 2;
    cache.f1linear = cache.gradf1x'*dir;
    cache.f1quad = cache.C1dir'*cache.QC1dir;
end

if prob.istheref2
    if prob.isthereC2
        if prob.isC2fun, cache.C2dir = prob.C2(dir);
        else cache.C2dir = prob.C2*dir; end
        ops.C2 = ops.C2 + 1;
    else
        cache.C2dir = dir;
    end
end

if prob.istherelin
    cache.lindir = prob.l'*dir;
end

cache.flagLineSearch = 1;
