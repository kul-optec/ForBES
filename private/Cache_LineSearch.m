function [cache, ops] = Cache_LineSearch(cache, dir)

if cache.flagLineSearch == 1
    if norm(dir-cache.dir) == 0
        ops = Ops_Init();
        return;
    end
end

if cache.flagGradStep == 0
    [cache, ops] = Cache_GradStep(cache, cache.gam);
else
    ops = Ops_Init();
end

prob = cache.prob;
cache.dir = dir;

if prob.istheref1
    if prob.isthereC1
        cache.C1dir = prob.C1*dir;
        cache.QC1dir = prob.Q(cache.C1dir);
        cache.C1tQC1dir = prob.C1'*cache.QC1dir;
        ops.C1 = ops.C1 + 2;
    else
        cache.C1dir = dir;
        cache.QC1dir = prob.Q(cache.C1dir);
        cache.C1tQC1dir = cache.QC1dir;
    end
    ops.gradf1 = ops.gradf1 + 1;
    cache.f1linear = cache.gradf1x(:)'*dir(:);
    cache.f1quad = cache.C1dir(:)'*cache.QC1dir(:);
end

if prob.istheref2
    if prob.isthereC2
        cache.C2dir = prob.C2*dir;
        ops.C2 = ops.C2 + 1;
    else
        cache.C2dir = dir;
    end
end

if prob.istherelin
    cache.lindir = prob.l(:)'*dir(:);
end

cache.flagLineSearch = 1;

end
