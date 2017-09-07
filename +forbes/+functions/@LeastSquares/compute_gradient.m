function [g, v] = compute_gradient(obj, x)
    res = obj.A*x - obj.b;
    g = obj.A'*res;
    v = 0.5*obj.lam*norm(res)^2;
end
