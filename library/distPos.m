%DISTPOS Distance from a box
%
%   DISTPOS(w, l) builds the function
%       
%       g(x) = sum(w_i*(x_i - max{l_i, x_i}))
%
%   Boundaries l_i can take the value -inf, in which case the corresponding
%   halfline is lower unbounded.
%
%   Weights w_i are assumed to be 1 if not provided. They can take the
%   value +inf, in which case the distance from the corresponding halfline
%   [l_i,+inf] becomes the indicator function.
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

function obj = distPos(weights, lb)
    if nargin < 1 || isempty(weights)
        weights = 1;
    end
    if nargin < 2 || isempty(lb)
        lb = 0;
    end
    if any(weights < 0)
        error('all weights must be nonnegative');
    end
    obj.makeprox = @() @(x, gam) call_distPos_prox(x, gam, lb, weights);
end

function [prox, val] = call_distPos_prox(x, gam, lb, weights)
    mu = gam*weights;
    prox = max(min(x+mu,lb),x);
    if nargout > 1
        finw = ~isinf(weights);
        val = sum(weights(finw).*max(lb(finw)-prox(finw),0));
    end
end
