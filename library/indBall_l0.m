%INDBALL_L0 Indicator function of the (nonconvex) L0 ball with given radius.
%
%   INDBALL_L0(N) builds the function
%
%       g(x) = 0    if nnz(x) <= N
%            = +inf otherwise
%
%   Argument N is required.

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

function obj = indBall_l0(N)
    obj.makeprox = @() @(x, gam) call_indBall_l0_prox(x, N);
end

function [prox, val] = call_indBall_l0_prox(x, N)
    prox = x;
    [~, I] = sort(prox, 'descend');
    prox(I(N+1:end)) = 0;
    val = 0;
end
