%LOGLOSS Allocates the log-logistic loss function.
%
%   LOGLOSS(mu) builds the log-logistic loss function
%
%       f(x) = (sum_i mu_i*log(1+exp(-x_i)))
%
%   If mu is a scalar then mu_i = mu for all i=1,...,n.
%   If mu is not provided then mu = 1.
%
% Copyright (C) 2015, Lorenzo Stella and Panagiotis Patrinos
%
% This file is part of ForBES.
%
% ForBES is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ForBES is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with ForBES. If not, see <http://www.gnu.org/licenses/>.

function obj = logLoss(mu)
    if nargin < 1, mu = 1; end
    if mu < 0, error('first argument mu must be nonnegative'); end
    obj.makef = @() @(x) call_logLoss_f(x, mu);
    obj.L = max(mu); % Lipschitz constant of the gradient of f
    obj.hasHessian = 1;
    obj.isConvex = 1;
end

function [val, grad, hess] = call_logLoss_f(x, mu)
  % value and gradient of f(x) = mu*sum(log(1+exp(-x)))
  emx = exp(-x);
  invpx = (1+emx);
  val = sum(mu.*log(invpx));
  if nargout >= 2
    px = 1./invpx;
    grad = mu.*(px-1);
    if nargout >= 3
      h = px.*(1-px);
      hess = diag(sparse(mu.*h));
    end
  end
end
