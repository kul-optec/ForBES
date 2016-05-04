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

function obj = dist2Neg(weights,ub)
    % Function value and gradient of (w/2)*dist^2(x,C) where C is the box [-infty,ub]
    if nargin < 1 || isempty(weights)
        weights = 1;
    end
    if nargin < 2 || isempty(ub)
        ub = 0;
    end
    if any(weights < 0)
        error('all weights must be nonnegative');
    end
    obj.makef = @() @(x) call_dist2Neg_f(x, ub, weights);
    obj.L = max(weights);
end

function [val, grad] = call_dist2Neg_f(x, ub, weights)
    proj = min(x,ub);
    diff = x - proj;
    grad = weights.*diff;
    val  = 0.5*sum(weights.*(diff.*diff));
end
