classdef ProblemComposite < handle
  properties
    x0, n, m1, m2
    f1, callf1, C1, d1, Q, q
    f2, callf2, C2, d2
    g, callg, D, mu
    lin
    istheref1, isthereC1
    istheref2, isthereC2
    isthereD
    istherelin
    Lf, exactLf, flagLipschitzConstant
    useHessian
  end
  methods
    function prob = ProblemComposite(f1, C1, d1, f2, C2, d2, g, D, lin, x0)
      prob.x0 = x0;
      prob.n = size(prob.x0);
      if isempty(f1)
        prob.istheref1 = false;
      else
        if ~isfield(f1, 'isQuadratic') || ~f1.isQuadratic
          error('function f1 must be quadratic');
        end
        if ~isfield(f1, 'makef')
          error('value/gradient of f1 is not defined (there is no makef)');
        end
        prob.istheref1 = true;
        prob.f1 = f1;
        prob.callf1 = f1.makef();
        if isempty(C1)
          prob.m1 = prob.n;
          prob.isthereC1 = false;
        else
          prob.C1 = C1;
          prob.m1 = [size(prob.C1, 1), 1];
          prob.isthereC1 = true;
        end
        [~, prob.q] = prob.callf1(zeros(prob.m1));
        prob.Q = @(x) prob.HessianQuadratic(x);
        if isempty(d1)
          prob.d1 = sparse(zeros(prob.m1));
        else
          prob.d1 = d1;
        end
      end
      if isempty(f2)
        prob.istheref2 = false;
      else
        if isfield(f2, 'isQuadratic') && f2.isQuadratic
          error('you should provide f2 as f1, since it is quadratic');
        end
        if ~isfield(f2, 'makef')
          error('value/gradient of f2 is not defined (there is no makef)');
        end
        prob.istheref2 = true;
        prob.f2 = f2;
        prob.callf2 = f2.makef();
        if isempty(C2)
          prob.m2 = prob.n;
          prob.isthereC2 = false;
        else
          prob.C2 = C2;
          prob.m2 = [size(prob.C2, 1), 1];
          prob.isthereC2 = true;
        end
        if isempty(d2)
          prob.d2 = zeros(prob.m2);
        else
          prob.d2 = d2;
        end
      end
      if isempty('g')
        error('missing g');
      end
      if ~isfield(g, 'makeprox')
        error('the prox for the term g you specified is not available');
      end
      prob.g = g;
      prob.callg = g.makeprox();
      if isempty(D)
        prob.isthereD = false;
      else
        [flag_tf, mu_D] = prob.isTightFrame(D);
        if ~flag_tf
          error('DD'' must be a positive multiple of the identity');
        end
        prob.mu = mu_D;
        prob.D = D;
        prob.isthereD = true;
      end
      if isempty(lin)
        prob.istherelin = false;
      else
        prob.istherelin = true;
        prob.lin = lin;
      end
      prob.useHessian = false;
    end
  end
end
