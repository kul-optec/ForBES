function [cache, ops] = Cache_LineSearch(cache, dir1)

if cache.flagGradStep == 0
    [cache, ops] = Cache_GradStep(cache, cache.gam);
else
    ops = Ops_Init();
end

prob = cache.prob;

if nargin > 1 && ~isempty(dir1)

    cache.dir1 = dir1;

    if prob.istheref1
        if prob.isthereC1
            cache.C1dir1 = prob.C1*dir1;
            cache.QC1dir1 = prob.Q(cache.C1dir1);
            cache.C1tQC1dir1 = prob.C1'*cache.QC1dir1;
            ops.C1 = ops.C1 + 2;
        else
            cache.C1dir1 = dir1;
            cache.QC1dir1 = prob.Q(cache.C1dir1);
            cache.C1tQC1dir1 = cache.QC1dir1;
        end
        ops.gradf1 = ops.gradf1 + 1;
        cache.f1linear = cache.gradf1x(:)'*dir1(:);
        cache.f1quad = cache.C1dir1(:)'*cache.QC1dir1(:);
    end

    if prob.istheref2
        if prob.isthereC2
            cache.C2dir1 = prob.C2*dir1;
            ops.C2 = ops.C2 + 1;
        else
            cache.C2dir1 = dir1;
        end
    end

    if prob.istherelin
        cache.lindir1 = prob.l(:)'*dir1(:);
    end

    cache.flagLineSearch1 = 1;

end

if nargin > 2 && ~isempty(dir2)

    cache.dir2 = dir2;

    if prob.istheref1
        if prob.isthereC1
            cache.C1dir2 = prob.C1*dir2;
            cache.QC1dir2 = prob.Q(cache.C1dir2);
            cache.C1tQC1dir2 = prob.C1'*cache.QC1dir2;
            ops.C1 = ops.C1 + 2;
        else
            cache.C1dir2 = dir2;
            cache.QC1dir2 = prob.Q(cache.C1dir2);
            cache.C1tQC1dir2 = cache.QC1dir2;
        end
        ops.gradf1 = ops.gradf1 + 1;
        cache.f1linear = cache.gradf1x(:)'*dir2(:);
        cache.f1quad = cache.C1dir2(:)'*cache.QC1dir2(:);
    end

    if prob.istheref2
        if prob.isthereC2
            cache.C2dir2 = prob.C2*dir2;
            ops.C2 = ops.C2 + 1;
        else
            cache.C2dir2 = dir2;
        end
    end

    if prob.istherelin
        cache.lindir2 = prob.l(:)'*dir2(:);
    end

    cache.flagLineSearch2 = 1;

end

end
