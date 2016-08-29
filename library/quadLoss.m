%QUADLOSS Allocates the squared norm function.
%
%   QUADLOSS(w, p) builds the function
%
%       f(x) = 0.5*sum_i w_i(x_i-p_i)^2
%
%   If w is a positive scalar then w_i = w (same for p). If omitted, w = 1
%   and p = 0.

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

function obj = quadLoss(w, p)
    if nargin < 1, w = 1; end
    if nargin < 2, p = 0; end
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
    elseif ismatrix(w)
        obj.makef = @() @(x) call_squaredWeightedDistance(x, w, p);
        if all(w > 0)
            obj.makefconj = @() @(x) call_squaredWeightedDistance_conj(x, 1./w, p);
        end
    end
end

function [v, g, H] = call_squaredWeightedDistance(x, w, p)
    res = x-p;
    g = w.*res;
    v = 0.5*(res(:)'*g(:));
    if nargout >= 3
        H = @(x) w.*x;
    end
end

function [v, g, H] = call_squaredWeightedDistance_conj(y, w, p)
    g = p + w.*y;
    v = 0.5*(y(:)'*(g(:) + p(:)));
    if nargout >= 3
        H = @(x) w.*x;
    end
end
