function [p, v] = compute_prox(obj, x, ~)
    if obj.method == 1 % using svds
        [U, S, V] = svds(x, obj.r, 'largest');
        p = U*(S*V');
    elseif obj.method == 2 % using lansvd
        [U, S, V] = lansvd(x, obj.r, 'L');
        p = U*(S*V');
    end
    v = 0;
end
