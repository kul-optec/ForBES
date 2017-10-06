% LASSO
%
%   LASSO(A, b, lam, opt) solves the problem
%
%     minimize (1/2)||Ax-b||^2 + lam*sum(abs(x))
%

function out = lasso(A, b, lam0, opt)

  % if options are not given
  if nargin < 4, opt = struct(); end

  % set some defualt options
  if ~isfield(opt, 'term'), opt.term = []; end
  term0 = opt.term;
  if ~isfield(opt, 'tol'), opt.tol = []; end
  tol0 = opt.tol;
  if ~isfield(opt, 'display'), opt.display = 0; end

  % compute Lipschitz constant
  [m, n] = size(A);
  eigsOpt.issym = 1;
  eigsOpt.tol = 1e-3;
  funHessian = @(x) A'*(A*x);
  Lf = eigs(funHessian, n, 1, 'LM', eigsOpt);
  opt.Lf = Lf;

  % to warm start or not to warm start?
  if ~isfield(opt, 'continuation') || isempty(opt.continuation), opt.continuation = 1; end
  if opt.continuation
    lam_max = norm(A'*b,'inf');
    lam = lam_max;
  else
    lam = lam0;
  end

  %
  f = quadLoss(1, zeros(m, 1));
  init = zeros(n, 1);

  for i_cont = 1:100

    % % this is the continuation scheme of SpaRSA
    % btilde = b-A*init;
    % lam = max(0.5*norm(A'*btilde,'inf'), lam0);

    % this is the simpler continuation scheme
    lam = max(0.5*lam, lam0);

    g = l1Norm(lam);
    if lam <= lam0
      opt.term = term0;
      opt.tol = tol0;
    else
      opt.term = [];
      opt.tol = 1e-3*lam;
    end
    out = forbes(f, g, init, {A, -b}, [], opt);
    if lam <= lam0
      break;
    end
    init = out.x;
  end

end
