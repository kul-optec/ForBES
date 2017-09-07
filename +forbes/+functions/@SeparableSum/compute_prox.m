function [p, v] = compute_prox(obj, x, gam)
    p = zeros(obj.dimsum(end), 1);
    v = 0;
    baseidx = 1;
    for i = 1:length(obj.idx)
        xcurr = reshape(x(baseidx:obj.dimsum(i)), obj.dims{i});
        [p1, v1] = obj.fs{obj.idx(i)}.compute_prox(xcurr, gam);
        p(baseidx:obj.dimsum(i)) = p1;
        v = v + v1;
        baseidx = obj.dimsum(i)+1;
    end
end
