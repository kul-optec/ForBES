%INDZERO Indicator function of the set {0}
%
%   INDZERO() builds the function
%
%       g(x) = 0    if x = 0
%            = +inf otherwise

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

function obj = indZero()
    obj.isConvex = 1;
    obj.makeprox = @() @(x, gam) call_indZero_proj(x);
end

function [prox, val] = call_indZero_proj(x)
    [n, m] = size(x);
    prox = zeros(n, m);
    val = 0;
end
