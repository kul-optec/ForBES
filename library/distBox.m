%DISTBOX Distance from a box
%
%   DISTBOX(l, u, w) builds the function
%       
%       g(x) = sum(w_i*(x_i - max{l_i, min{u_i, x_i}}))
%
%   Boundaries l_i and u_i can take the value -inf and +inf respectively,
%   in which case the corresponding segment is lower or upper unbounded.
%
%   Weights w_i are assumed to be 1 if not provided. They can take the
%   value +inf, in which case the distance from the corresponding segment
%   [l_i,u_i] becomes the indicator function.
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

function obj = distBox(lb, ub, weights)
    if nargin < 3 || isempty(weights)
        weights = 1;
    end
    if nargin < 2 || isempty(ub)
        ub = +inf;
    end
    if nargin < 1 || isempty(lb)
        lb = -inf;
    end
    if any(weights < 0)
        error('all weights must be nonnegative');
    end
    obj.makeprox = @() @(x, gam) call_distBox_prox(x, gam, lb, ub, weights);
end

function [prox, val] = call_distBox_prox(x, gam, lb, ub, weights)
    mu = gam*weights;
    prox = max(x-ub-mu, 0) - max(lb-x-mu, 0) + min(max(x, lb), ub);
    if nargout > 1
        finw = ~isinf(weights);    
        val = sum(weights(finw).*abs(prox(finw)-min(max(prox(finw),lb(finw)),ub(finw))));
    end
end
