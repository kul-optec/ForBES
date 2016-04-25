%INDBIN Allocates the indicator function of a binary variable.
%
%   INDBIN(v0, v1) builds indicator function of the set {v0, v1}^n.
%   By default, v0 = 0 and v1 = 1.
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

function obj = indBin(v0, v1)
    if nargin < 1, v0 = 0; end
    if nargin < 2, v1 = 1; end
    obj.makeprox = @() @(x, gam) call_indBin_prox(x, v0, v1);
    obj.isConvex = 1;
end

% TODO: implement faster solution for particular cases
% like {-1, 1} (using sign) or {0, 1} (using >= 0.5)

function [prox, g] = call_indBin_prox(x, v0, v1)
    mid = (v1+v0)/2;
    ind = (x >= mid);
    prox(ind, 1) = v1;
    prox(~ind, 1) = v0;
    g = 0;
end
