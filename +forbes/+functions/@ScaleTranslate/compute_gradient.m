function [grad, val] = gradient(obj, x)
    [grad1, val] = obj.f.compute_gradient(obj.a*x + obj.b);
    grad = obj.a*grad1;
end
