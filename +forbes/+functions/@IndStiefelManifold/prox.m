function [p, v] = prox(obj, x, gam)
    if size(x, 1) < size(x, 2)
        error('argument must be a tall matrix');
    end
    [U, ~, V] = svd(x, 'econ');
    p = U*V';
    v = 0;
end
