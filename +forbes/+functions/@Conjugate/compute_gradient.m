function [g, v] = compute_gradient(obj, x)
    [g, v] = obj.f.compute_gradient_conjugate(x);
end
