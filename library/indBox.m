%INDBOX Indicator function of a box.
%
%   INDBOX(l, u) builds the function
%       
%       g(x) = 0    if l_i <= x_i <= u_i for all i
%            = +inf otherwise
%
%   Arguments l and u can be either vectors of the same size of x
%   or scalars. If any l_i (u_i) has value -inf (+inf) then the
%   corresponding segment [l_i, u_i] is lower (upper) unbounded.
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

function obj = indBox(lower, upper)
    obj.makeprox = @() @(x, gam) call_indBox_prox(x, lower, upper);
end

function [prox, val] = call_indBox_prox(x, lower, upper)
    prox = min(upper, max(lower, x));
    val = 0;
end
