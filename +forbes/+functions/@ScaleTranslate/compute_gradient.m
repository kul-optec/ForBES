function [grad, val] = compute_gradient(obj, x)
    [grad1, val] = obj.f.compute_gradient(obj.a*x + obj.b);
    grad = obj.a*grad1;
end
