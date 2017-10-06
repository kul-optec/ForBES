function [grad, val, hess] = compute_gradient(obj, x)
    emx = exp(-x);
    invpx = (1+emx);
    val = sum(log(invpx))*obj.mu;
    if nargout >= 2
        px = 1./invpx;
        grad = (px-1)*obj.mu;
        if nargout >= 3
            h = px.*(1-px);
            hess = obj.mu*diag(sparse(h));
        end
    end
end
