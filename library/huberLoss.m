%HUBERLOSS Allocates the Huber loss function.
%
%   HUBERLOSS(del) builds the function
%       
%       f(x) = sum_i l(x_i)
%
%   where
%
%       l(x_i) = 0.5/del*x_i^2      if |x_i| <= del
%                |x_i|-0.5*del      otherwise 
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
    %
    % Only f available for this function
    %
    obj.makef = @() @(x) call_huberLoss_f(x, del);
    obj.L = 1/del;
end

function [val, grad] = call_huberLoss_f(x, del)
    absx = abs(x);
    small = absx <= del;
    large = ~small;
    sqx = (0.5/del)*(x(small).^2);
    linx = absx(large)-0.5*del;
    val = sum(sqx)+sum(linx);
    if nargout >= 2
        grad = zeros(length(x),1);
        grad(small) = x(small)/del;
        grad(large) = sign(x(large));
    end
end
