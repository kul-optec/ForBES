function [g, v] = gradient(obj, x)
    [g, v] = obj.f.gradient_conjugate(x);
end
