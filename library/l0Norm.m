%L0NORM Allocates the L0 norm function.
%
%   L0NORM(mu) builds the function
%
%       g(x) = mu*||x||_0 = mu*nnz(x)

% Copyright (C) 2015-2016, Lorenzo Stella and Panagiotis Patrinos
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
    if nargin < 1
        mu = 1;
    end
    obj.makeprox = @() @(x, gam) call_l0Norm_prox(x, gam, mu);
end

function [prox, g] = call_l0Norm_prox(x, gam, mu)
    over = abs(x) > sqrt(2*gam*mu);
    prox = x.*over;
    if nargout >= 2
        g = mu*nnz(prox);
    end
end
