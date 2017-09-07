function [grad, val] = compute_gradient(obj, x)
    absx = abs(x);
    small = absx <= obj.del;
    large = ~small;
    sqx = (0.5/obj.del)*(x(small).^2);
    linx = absx(large)-0.5*obj.del;
    val = sum(sqx)+sum(linx);
    if nargout >= 2
        grad = zeros(length(x),1);
        grad(small) = x(small)/obj.del;
        grad(large) = sign(x(large));
    end
end
