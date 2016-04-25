%L1NORM Allocates the L1 norm function.
%
%   L1NORM(mu) builds the function
%       
%       g(x) = sum_i mu_i*abs(x_i)
%
%   If mu is a scalar then mu_i = mu for all i = 1,...,n.
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

function obj = l1Norm(mu)
    if nargin < 1, mu = 1; end
    if mu < 0, error('first argument mu must be nonnegative'); end
    obj.makeprox = @() @(x, gam) call_l1Norm_prox(x, gam, mu);
    obj.isConvex = 1;
end

function [prox, g] = call_l1Norm_prox(x, gam, mu)
    uz = max(0, abs(x)-gam*mu);
    prox = sign(x).*uz;
    if nargout >= 2
        g = sum(mu.*uz);
    end
end
