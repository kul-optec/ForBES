function [L, exact] = EstimateLipschitz(prob)

if prob.flagLipschitzConstant
  L = prob.Lf;
  exact = prob.exactLf;
  return;
end


% if prob.istheref1 && ~prob.istheref2 && prob.f1.isConvex
%   % if the smooth term is purely quadratic
%   % then compute L 'exactly'
%   eigsOpt.issym = 1;
%   eigsOpt.tol = 1e-3;
%   
%   if prob.isthereC1
%     funHessian = @(x) vec(prob.C1'*(prob.Q(prob.C1*reshape(x, prob.n))));
%     L          = eigs(funHessian, prod(prob.n), 1, 'LM', eigsOpt);
%   else
%     funHessian = @(x) vec(prob.Q(reshape(x, prob.n)));
%     L          = eigs(funHessian, prod(prob.n), 1, 'LM', eigsOpt);
%   end
%   
%   exact        = true;
%   prob.Lf      = L;
%   prob.exactLf = exact;
%   prob.flagLipschitzConstant = true;
%   return;
% end

delta = max(1e-12, prob.x0*1e-6);

if prob.istheref1
  if prob.isthereC1
    C1x0 = prob.C1*prob.x0;
    C1x1 = prob.C1*(prob.x0+delta);
    [~, gradf1C1x0] = prob.callf1(C1x0);
    [~, gradf1C1x1] = prob.callf1(C1x1);
    grad1x0 = prob.C1'*gradf1C1x0;
    grad1x1 = prob.C1'*gradf1C1x1;
  else
    [~, grad1x0] = prob.callf1(prob.x0);
    [~, grad1x1] = prob.callf1(prob.x0+delta);
  end
else
  grad1x0 = 0;
  grad1x1 = 0;
end

if prob.istheref2
  if prob.isthereC2
    C2x0 = prob.C2*prob.x0;
    C2x1 = prob.C2*(prob.x0+delta);
    [~, gradf2C2x0] = prob.callf2(C2x0);
    [~, gradf2C2x1] = prob.callf2(C2x1);
    grad2x0 = prob.C2'*gradf2C2x0;
    grad2x1 = prob.C2'*gradf2C2x1;
  else
    [~, grad2x0] = prob.callf2(prob.x0);
    [~, grad2x1] = prob.callf2(prob.x0+delta);
  end
else
  grad2x0 = 0;
  grad2x1 = 0;
end

exact   = true;
L       = norm(grad1x1 + grad2x1 - grad1x0 - grad2x0)/norm(delta);
prob.Lf = L;
prob.exactLf = exact;
prob.flagLipschitzConstant = true;
