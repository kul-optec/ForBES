%L0NORM Allocates the L0 norm function.
%
%   L0NORM(mu) builds the weighted L0 norm function
%       
%       g(x) = sum_i mu_i*(x_i ~= 0)
%
%   If mu is a scalar then mu_i = mu for all i=1,...,n and
%
%       g(x) = mu*||x||_0 = mu*nnz(x)
%
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

function obj = l0Norm(mu)
    if nargin < 1, mu = 1; end
    obj.makeprox = @() @(x, gam) call_l0Norm_prox(x, gam, mu);
    obj.isConvex = 0;
end

function [prox, g] = call_l0Norm_prox(x, gam, mu)
    over = abs(x) > sqrt(2*gam*mu); 
    prox = x.*over;
    if nargout >= 2
        g = mu*nnz(prox);
    end
end
