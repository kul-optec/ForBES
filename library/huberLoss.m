%HUBERLOSS Allocates the Huber loss function.
%
%   HUBERLOSS(del) builds the function
%       
%       f(x) = sum_i l_i(x_i)
%
%   where
%
%       l_i(x_i) = (0.5/del_i)*x_i^2      if |x_i| <= del_i
%                  |x_i|-0.5*del_i        otherwise 
%
%   If del is a scalar then del_i = del for i = 1,...,n.
%   If del is not provided then del = 1.
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

function obj = huberLoss(del)
    if nargin < 1 || isempty(del), del = 1; end
    if any(del <= 0)
        error('first argument del must be positive');
    end
    obj.makef = @() @(x) call_huberLoss_f(x, del);
    obj.L = 1/min(del);
    obj.isConvex = 1;
end

function [val, grad] = call_huberLoss_f(x, del)
    absx = abs(x);
    small = absx <= del;
    large = ~small;
    sqx = (0.5./del).*(x(small).^2);
    linx = absx(large)-0.5*del;
    val = sum(sqx)+sum(linx);
    if nargout >= 2
        grad = zeros(length(x),1);
        grad(small) = x(small)./del(small);
        grad(large) = sign(x(large));
    end
end
