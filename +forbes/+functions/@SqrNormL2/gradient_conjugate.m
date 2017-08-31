function [g, v] = gradient_conjugate(obj, x)
    if any(obj.w == 0)
        error('the conjugate is not smooth');
    end
    g = x ./ obj.w;
    v = 0.5*(x(:)'*g(:));
    if nargout >= 3
        H = @(x) x ./ obj.w;
    end
end
