function [p, v] = prox(obj, x, ~)
    [m_x, n_x] = size(x);
    if obj.method == 1 % using svds
        [U, S, V] = svds(reshape(x, obj.m, obj.n), obj.r, 'largest');
        p = reshape(U*(S*V'), m_x, n_x);
    elseif obj.method == 2 % using lansvd
        [U, S, V] = lansvd(reshape(x, obj.m, obj.n), obj.r, 'L');
        p = reshape(U*(S*V'), m_x, n_x);
    end
    v = 0;
end
