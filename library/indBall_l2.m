%INDBALL_L2 Indicator function of the L2 ball with given center and radius.
%
%   INDBALL_L2(rho, c) builds the function
%       
%       g(x) = 0    if ||x-c|| <= rho
%            = +inf otherwise
%
%   If c is not provided, c = 0. If also rho is not provided, rho = 1.
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

function obj = indBall_l2(rho, c)
    if nargin < 1 || isempty(rho)
        rho = 1;
    end
    if nargin < 2 || isempty(c)
        c = 0;
    end
    obj.makeprox = @() @(x, gam) call_indBall_l2_prox(x, rho, c);
end

function [prox, val] = call_indBall_l2_prox(x, rho, c)
    xc = x - c;
    nxc = norm(xc);
    if nxc <= rho
        prox = x;
    else
        prox = c + (rho/nxc)*xc;
    end
    val = 0;
end
