function [g, v, Q] = compute_gradient(obj, x)
    g = obj.Q*x+obj.q;
    v = 0.5*(g+obj.q)'*x;
    Q = obj.Q;
end
