function Set_Directions(cache, dir1, dir2)

if cache.flagGradStep == 0
    cache.Get_GradStep();
end

prob = cache.prob;

if nargin > 1 && ~isempty(dir1)

    cache.dir1 = dir1;

    if prob.istheref1
        if prob.isthereC1
            cache.C1dir1 = prob.C1*dir1;
            cache.QC1dir1 = prob.Q(cache.C1dir1);
            cache.C1tQC1dir1 = prob.C1'*cache.QC1dir1;
            cache.ops.addC1();
            cache.ops.addC1(); % yes, twice
        else
            cache.C1dir1 = dir1;
            cache.QC1dir1 = prob.Q(cache.C1dir1);
            cache.C1tQC1dir1 = cache.QC1dir1;
        end
        cache.ops.addgradf1();
        cache.f1linear1 = cache.gradf1x(:)'*dir1(:);
        cache.f1quad1 = cache.C1dir1(:)'*cache.QC1dir1(:);
    end

    if prob.istheref2
        if prob.isthereC2
            cache.C2dir1 = prob.C2*dir1;
            cache.ops.addC2();
        else
            cache.C2dir1 = dir1;
        end
    end

    if prob.istherelin
        cache.lindir1 = prob.lin(:)'*dir1(:);
    end

    cache.flagLineSearch1 = true;

end

if nargin > 2 && ~isempty(dir2)

    cache.dir2 = dir2;

    if prob.istheref1
        if prob.isthereC1
            cache.C1dir2 = prob.C1*dir2;
            cache.QC1dir2 = prob.Q(cache.C1dir2);
            cache.C1tQC1dir2 = prob.C1'*cache.QC1dir2;
            cache.ops.addC1();
            cache.ops.addC1(); % yes, twice
        else
            cache.C1dir2 = dir2;
            cache.QC1dir2 = prob.Q(cache.C1dir2);
            cache.C1tQC1dir2 = cache.QC1dir2;
        end
        cache.ops.addgradf1();
        cache.f1linear2 = cache.gradf1x(:)'*dir2(:);
        cache.f1quad2 = cache.C1dir2(:)'*cache.QC1dir2(:);
        cache.f1cross = cache.QC1dir1(:)'*cache.C1dir2(:);
    end

    if prob.istheref2
        if prob.isthereC2
            cache.C2dir2 = prob.C2*dir2;
            cache.ops.addC2();
        else
            cache.C2dir2 = dir2;
        end
    end

    if prob.istherelin
        cache.lindir2 = prob.lin(:)'*dir2(:);
    end

    cache.flagLineSearch2 = true;

end
