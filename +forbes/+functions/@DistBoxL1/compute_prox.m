function [p, v] = compute_prox(obj, x, gam)
    mu = gam * obj.w;
    p = max(x-obj.ub-mu, 0) - max(obj.lb-x-mu, 0) + min(max(x, obj.lb), obj.ub);
    if nargout > 1
        proj = max(obj.lb, min(obj.ub, p));
        if isscalar(obj.w)
            v = sum(obj.w*abs(p-proj));
        else
            finw = ~isinf(obj.w);
            v = sum(obj.w(finw).*abs(p(finw)-proj(finw)));
        end
    end
end
