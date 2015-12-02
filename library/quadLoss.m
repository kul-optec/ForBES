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
    if isscalar(w)
        obj.Q = w;
        obj.q = -w*p;
    elseif isvector(w)
        n = length(w);
        obj.Q = spdiags(w,0,n,n);
        obj.q = -w.*p;
    end
    obj.isQuadratic = 1;
    obj.isConjQuadratic = 1;
    if all(w > 0)
        obj.makefconj = @() @(x) call_squaredWeightedDistance_conj(x, w, p);
    end
end

function [val, grad] = call_squaredWeightedDistance_conj(y, w, p)
    grad = p + (y./w);
    val = 0.5*(y'*(grad + p));
end
