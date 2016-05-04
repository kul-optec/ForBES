function obj = stackOp(objs, idx)
    l = length(objs);
    if nargin < 2
        idx = 1:l;
    end
    n = objs{1}.n;
    msum = zeros(1, length(idx));
    msum(1) = objs{1}.m;
    for i = 2:length(idx)
        if objs{idx(i)}.n ~= n, error('all operators must be defined over the same space'); end
        msum(i) = msum(i-1)+objs{idx(i)}.m;
    end
    callops = {};
    calladjs = {};
    for i = 1:l
        callops{i} = objs{i}.makeop();
        calladjs{i} = objs{i}.makeadj();
    end
    obj.m = msum(end);
    obj.n = n;
    obj.makeop = @() @(x) call_stackOps(x, callops, idx);
    obj.makeadj = @() @(y) call_sumOps(y, calladjs, idx, msum);
end

function y = call_stackOps(x, ops, idx)
    y = [];
    for i = 1:length(idx)
        y = [y; ops{idx(i)}(x)];
    end
end

function x = call_sumOps(y, ops, idx, msum)
    x = 0;
    baseidx = 1;
    for i = 1:length(idx)
        x = x + ops{idx(i)}(y(baseidx:msum(i)));
        baseidx = msum(i)+1;
    end
end
