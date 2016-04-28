%L2NORMSUM Allocates the sum-of-L2-norm function.
%
%   L2NORMSUM(m, mu) builds the function
%       
%       g(x) = mu*sum(||x_i||_2)
%
%   where x_i are blocks of size m of vector x. If mu is not provided it is
%   assumed mu = 1. If also m is not provided, then m = 1 and the function
%   is the L1 norm.
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

function obj = l2NormSum(m, mu)
    if nargin < 2
        mu = 1;
        if nargin < 1
            m = 1;
        end
    end
    obj.makeprox = @() @(x, gam) call_l2NormSum_prox(x, gam, m, mu);
end

function [z, v] = call_l2NormSum_prox(x, gam, m, mu)
    n = length(x);
    nb = n/m;
    x_resh = reshape(x, m, nb);
    modx = sqrt(sum(x_resh.*x_resh, 1));
    newmodx = max(0, modx-gam*mu);
    scal_block = newmodx./max(gam*mu, modx);
    scal_resh = repmat(scal_block, m, 1);
    scal = reshape(scal_resh, n, 1);
    z = scal.*x;
    if nargout >= 2
        v = mu*sum(newmodx);
    end
end
