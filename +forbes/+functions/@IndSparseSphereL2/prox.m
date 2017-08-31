function [p, v] = prox(obj, x, ~)
    p = x;
    [~, I] = sort(abs(p), 'descend');
    p(I(obj.N+1:end)) = 0;
    normp = norm(p);
    p = obj.R*p/normp;
    v = 0;
end
