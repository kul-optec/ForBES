%INDNEG Indicator function of the negative orthant.
%
%   INDNEG(u) builds the function
%       
%       g(x) = 0    if x_i <= u_i for all i
%            = +inf otherwise
%
%   Argument u is either a scalar or a vector of the same size of x. If
%   argument u is not given, then u = 0. If any u_i is +inf then the
%   corresponding halfline [-inf, u_i] is unbounded.
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

function obj = indNeg(ub)
    if nargin < 1 || isempty(ub)
        ub = 0;
    end
    obj.makeprox = @() @(x, gam) call_indNeg_prox(x, ub);
end

function [prox, val] = call_indNeg_prox(x, ub)
    prox = min(ub, x);
    val = 0;
end
