function [p, v] = compute_prox(obj, x, gam)
    p = x;
    [~, I] = sort(abs(p), 'descend');
    p(I(obj.N+1:end)) = 0;
    v = 0;
end