function Hx = Util_HessianQuadratic(x, f, grad0)

[~, grad] = f(x);
Hx = grad-grad0;

end
