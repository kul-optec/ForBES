function obj = sumOp(objs, idx)
    l = length(objs);
    if nargin < 2
        idx = 1:l;
    end
    nsum = zeros(1, length(idx));
    nsum(1) = objs{1}.n;
    m = objs{1}.m;
    for i = 2:length(idx)
        if objs{idx(i)}.m ~= m, error('all operators must have image in the same space'); end
        nsum(i) = nsum(i-1) + objs{idx(i)}.n;
    end
    callops = {};
    calladjs = {};
    for i = 1:l
        callops{i} = objs{i}.makeop();
        calladjs{i} = objs{i}.makeadj();
    end
    obj.m = m;
    obj.n = nsum(end);
    obj.makeop = @() @(x) call_sumOps(x, callops, idx, nsum);
    obj.makeadj = @() @(y) call_stackOps(y, calladjs, idx);
end

function x = call_stackOps(y, ops, idx)
    x = [];
    baseidx = 0;
    for i = 1:length(idx)
        x = [x; ops{idx(i)}(y)];
    end
end

function y = call_sumOps(x, ops, idx, nsum)
    y = 0;
    baseidx = 1;
    for i = 1:length(idx)
        y = y + ops{idx(i)}(x(baseidx:nsum(i)));
        baseidx = nsum(i)+1;
    end
end
