function Hx = HessianQuadratic(prob, x)

[~, grad] = prob.callf1(x);
Hx = grad-prob.q;
