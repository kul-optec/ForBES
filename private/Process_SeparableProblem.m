function [prob, dualprob] = Process_SeparableProblem(prob, opt)

if isfield(prob, 'b')
    if norm(prob.b) > 0,
        dualprob.l = prob.b;
        dualprob.istherelin = true;
    else
        dualprob.istherelin = false;
    end
    n = size(prob.b);
    if isfield(prob, 'y0'), dualprob.x0 = prob.y0;
    else dualprob.x0 = zeros(n); end
    dualprob.n = n;
else
    error('you must specify the right hand side b of the equality constraint');
end
if ~any(isfield(prob, {'f1', 'f2'}))
    error('missing f1 and f2');
end
if isfield(prob, 'f1')
    if ~isfield(prob.f1, 'isConjQuadratic') || ~prob.f1.isConjQuadratic
        error('the conjugate function of f1 must be quadratic');
    end
    if ~isfield(prob.f1, 'makefconj'), error('conjugate function of f1 is not defined'); end
    dualprob.callf1 = prob.f1.makefconj();
    dualprob.istheref1 = true;
    if isfield(prob, 'A1')
        dualprob.isthereC1 = true;
        dualprob.isC1fun = false;
        dualprob.C1 = -prob.A1';
        dualprob.m1 = [size(dualprob.C1, 1), 1];
        dualprob.d1 = zeros(dualprob.m1);
    else
        error('you must specify matrix A1 in the constraint');
    end
    [~, q] = dualprob.callf1(zeros(dualprob.m1));
    dualprob.Q = @(x) Util_HessianQuadratic(x, dualprob.callf1, q);
else
    dualprob.istheref1 = false;
end
if isfield(prob, 'f2')
    if isfield(prob.f2, 'isConjQuadratic') && prob.f2.isConjQuadratic
        error('you should provide f2 as f1, since its conjugate is quadratic');
    end
    dualprob.istheref2 = true;
    if ~isfield(prob.f2, 'makefconj'), error('conjugate function of f2 is not defined'); end
    dualprob.callf2 = prob.f2.makefconj();
    if isfield(prob, 'A2')
        dualprob.isthereC2 = true;
        dualprob.isC2fun = false;
        dualprob.C2 = -prob.A2';
        dualprob.m2 = [size(dualprob.C2, 1), 1];
        dualprob.d2 = zeros(dualprob.m2);
    else
        error('yout must specify matrix A2 in the constraint');
    end
    if isfield(prob.f2, 'mu'), dualprob.Lf2 = 1/prob.f2.mu; end
    dualprob.useHessian = 0;
else
    dualprob.istheref2 = false;
end
if isfield(prob, 'B')
    if isa(prob.B, 'function_handle'), error('you must specify matrix B as a matrix in the constraint'); end
else
    error('you must specify matrix B in the constraint');
end
mus = sum((prob.B).*(prob.B), 1);
if (max(mus)-min(mus))/max(mus) > 10*eps, error('B''B must be a positive multiple of the identity'); end
prob.muB = mus(1);
if ~isfield(prob, 'g'), error('you must specify term g'); end
if ~isfield(prob.g, 'makeprox'), error('the prox for the term g you specified is not available'); end
prob.callg = prob.g.makeprox();
dualprob.callg = make_prox_conj(prob.callg, prob.B, prob.muB);
if isfield(opt, 'Lf')
    dualprob.Lf = opt.Lf;
else
    dualprob.Lf = Process_LipschitzConstant(dualprob);
end
dualprob.muf = 0;
dualprob.processed = true;

end

function op = make_prox_conj(proxg, B, mu)

op = @(y, gam) call_prox_conj(y, gam, proxg, B, mu);

end

function [proxpoint, proxval] = call_prox_conj(y, gam, prox, B, mu)
% Compute prox_{gamma,h}, where h(y) = g*(-B'y), from prox_{gamma,g} and B.
% Uses Bauschke, Combettes, Prop. 23.32, and Moreau identity.
% Then uses conjugate subgradient theorem to compute the function value.

mugam = mu*gam;
[z, v] = prox(-(B'*y)/mugam, 1/mugam);
Bz = B*z;
proxpoint = y+gam*Bz;
proxval = -proxpoint'*Bz - v;

end
