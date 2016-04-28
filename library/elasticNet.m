%ELASTICNET Allocates the elastic net regularization function.
%
%   ELASTICNET(mu, lam) builds the function
%       
%       g(x) = mu*||x||_1 + (lam/2)*||x||^2
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
    if nargin < 2
        lam = 1;
        if nargin < 1, mu = 1; end
    end
    obj.makeprox = @() @(x, gam) call_elasticNet_prox(x, gam, mu, lam);
end

function [prox, g] = call_elasticNet_prox(x, gam, mu, lam)
    uz = max(0, abs(x)-gam*mu)/(1+lam*gam);
    prox = sign(x).*uz;
    if nargout >= 2
        g = mu*sum(uz)+(0.5*lam)*(uz'*uz);
    end
end
