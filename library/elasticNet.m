%ELASTICNET Allocates the elastic net regularization function.
%
%   ELASTICNET(mu, lam) builds the function
%       
%       g(x) = sum_i mu_i*abs(x_i) + (lam_i/2)*x_i^2
%
%   If mu (or lam) is a scalar then mu_i = mu (or lam_i = lam) for all
%   i = 1,...,n. If mu is not provided or is empty then mu = 1 (the same
%   holds for lam).
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

function obj = elasticNet(mu, lam)
    if nargin < 1 || isempty(mu), mu = 1; end
    if nargin < 2 || isempty(lam), lam = 1; end
    if any(mu < 0), error('first argument mu must be nonnegative'); end
    if any(lam < 0), error('second argument lam must be nonnegative'); end
    obj.makeprox = @() @(x, gam) call_elasticNet_prox(x, gam, mu, lam);
    obj.isConvex = 1;
end

function [prox, g] = call_elasticNet_prox(x, gam, mu, lam)
    uz = max(0, abs(x)-gam*mu)./(1+lam*gam);
    prox = sign(x).*uz;
    if nargout >= 2
        g = sum(mu.*uz)+0.5*(uz'*(lam.*uz));
    end
end
