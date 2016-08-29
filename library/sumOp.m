function obj = sumOp(objs, idx)
    l = length(objs);
    if nargin < 2, idx = 1:l; end
    
    nsum = zeros(1, length(idx));
    dims = [zeros(length(idx), 1), ones(length(idx), 1)];
    nsum(1) = prod(objs{1}.n);
    dims(1,1) = objs{1}.n(1);
    if size(objs{1}.n, 2) > 1, dims(1,2) = objs{1}.n(2); end
    m = objs{1}.m;
    for i = 2:length(idx)
        if any(objs{idx(i)}.m ~= m), error('all operators must have image in the same space'); end
        dims(i,1) = objs{i}.n(1);
        if size(objs{i}.n, 2) > 1, dims(i,2) = objs{i}.n(2); end
        nsum(i) = nsum(i-1) + prod(dims(i,:));
    end
    
    callops = {};
    calladjs = {};
    for i = 1:l
        callops{i} = objs{i}.makeop();
        calladjs{i} = objs{i}.makeadj();
    end
    
    obj.m = m;
    obj.n = [nsum(end), 1];
    obj.makeop = @() @(x) call_sumOps(x, callops, idx, nsum, dims);
    obj.makeadj = @() @(y) call_stackOps(y, calladjs, idx);
end

function y = call_sumOps(x, ops, idx, nsum, dims)
    y = 0;
    baseidx = 1;
    for i = 1:length(idx)
        y = y + ops{idx(i)}(reshape(x(baseidx:nsum(i)), dims(i,1), dims(i,2)));
        baseidx = nsum(i)+1;
    end
end

function x = call_stackOps(y, ops, idx)
    x = [];
    for i = 1:length(idx)
        z = ops{idx(i)}(y);
        x = [x; z(:)];
    end
end
