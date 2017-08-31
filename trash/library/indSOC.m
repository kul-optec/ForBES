%INDSOC Indicator function of the positive orthant.
%
%   INDSOC(b) builds the indicator of the translated second-order cone, i.e.
%
%       g(x) = 0    if norm(z(2:end)) <= z(1), where z = x-b
%            = +inf otherwise
%
%   If argument b is not given or empty, then b = 0.

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

function obj = indSOC(b)
    if nargin < 1 || isempty(b), b = 0; end
    obj.isConvex = 1;
    obj.makeprox = @() @(x, gam) call_indSOC_prox(x, b);
end

function [prox, val] = call_indSOC_prox(x, b)
    z = x-b;
    
    if isempty(z)
        z=[];
        return;
    elseif length(z)==1
        z = max(z,0);
        return;
    end

    v1 = z(1);
    v2 = z(2:end);
    normv2 = norm(v2);

    if v1 <= -normv2
        z = zeros(length(z), 1);
    elseif v1 >= normv2
        z = z;
    else
        a = (v1+normv2)/2;
        z(1) = a;
        z(2:end) = a*(z(2:end)/normv2);
    end
    
    prox = z+b;
    val = 0;
end
