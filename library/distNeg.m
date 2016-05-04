%DISTNEG Distance from a box
%
%   DISTNEG(w, u) builds the function
%       
%       g(x) = sum(w_i*(x_i - min{u_i, x_i}))
%
%   Boundaries u_i can take the value +inf, in which case the corresponding
%   halfline is upper unbounded.
%
%   Weights w_i are assumed to be 1 if not provided. They can take the
%   value +inf, in which case the distance from the corresponding halfline
%   [-inf,u_i] becomes the indicator function.
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

function obj = distNeg(weights, ub)
    if nargin < 1 || isempty(weights)
        weights = 1;
    end
    if nargin < 2 || isempty(ub)
        ub = 0;
    end
    if any(weights < 0)
        error('all weights must be nonnegative');
    end
    obj.makeprox = @() @(x, gam) call_distNeg_prox(x, gam, ub, weights);
end

function [prox, val] = call_distNeg_prox(x, gam, ub, weights)
    mu = gam*weights;
    prox = min(max(x-mu,ub),x);
    if nargout > 1
        finw = ~isinf(weights);
        val = sum(weights(finw).*max(prox(finw)-ub(finw),0));
    end
end
