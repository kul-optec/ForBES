function [g, v] = compute_gradient(obj, x)
    g = zeros(obj.dimsum(end), 1);
    v = 0;
    baseidx = 1;
    for i = 1:length(obj.idx)
        xcurr = reshape(x(baseidx:obj.dimsum(i)), obj.dims{i});
        [g1, v1] = obj.fs{obj.idx(i)}.compute_gradient(xcurr);
        g(baseidx:obj.dimsum(i)) = g1;
        v = v + v1;
        baseidx = obj.dimsum(i)+1;
    end
end
