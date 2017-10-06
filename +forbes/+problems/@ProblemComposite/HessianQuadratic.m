function Hx = HessianQuadratic(prob, x)

[grad, ~] = prob.f1.gradient(x);
Hx = grad-prob.q;
