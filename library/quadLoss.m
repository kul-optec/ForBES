%QUADLOSS Allocates the squared norm function.
%
%   QUADLOSS(w, p) builds the function
%       
%       f(x) = 0.5*sum_i w_i(x_i-p_i)^2
%   
%   Both arguments are required. If w is a positive scalar then w_i = w.
%   Length of p must instead be compliant with the dimension of the domain
%   of the function.
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

function obj = quadLoss(w, p)
    if any(w < 0)
        error('first argument should be nonnegative');
    end
    obj.isConvex = 1;
    obj.isQuadratic = 1;
    obj.isConjQuadratic = 1;
    obj.hasHessian = 1;
    obj.L = max(w);
    if isscalar(w)
        obj.makef = @() @(x) call_squaredWeightedDistance(x, w, p);
        if w > 0
            obj.makefconj = @() @(x) call_squaredWeightedDistance_conj(x, 1/w, p);
        end
    elseif isvector(w)
        n = length(w);
        obj.makef = @() @(x) call_squaredWeightedDistance(x, spdiags(w,0,n,n), p);
        if all(w > 0)
            obj.makefconj = @() @(x) call_squaredWeightedDistance_conj(x, spdiags(1./w,0,n,n), p);
        end
    end
end

function [v, g, W] = call_squaredWeightedDistance(x, W, p)
    res = x-p;
    g = W*res;
    v = 0.5*(res'*g);
end

function [v, g, W] = call_squaredWeightedDistance_conj(y, W, p)
    g = p + W*y;
    v = 0.5*(y'*(g + p));
end
