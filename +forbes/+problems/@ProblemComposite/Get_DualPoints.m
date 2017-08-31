function [x1, x2, z] = Get_DualPoints(prob, y, gam)

Ax = 0;
if prob.istheref1
    [~, x1] = prob.callf1(prob.C1*y);
    Ax = Ax - prob.C1'*x1;
else
    x1 = [];
end
if prob.istheref2
    [~, x2] = prob.callf2(prob.C2*y);
    Ax = Ax - prob.C2'*x2;
else
    x2 = [];
end
w = -gam*prob.D*(prob.lin - Ax - y/gam);
u = -prob.D'*(prob.callg(w, prob.mu*gam) - w)/(prob.mu*gam);
z = (prob.D*u)/prob.mu;
%w + prob.D'*(proxg(prob.D*w, prob.mu/gam) - prob.D*w)/prob.mu;
%-prob.D'*prob.callg((prob.D*w)/(prob.mu*gam), 1/(prob.mu*gam));
