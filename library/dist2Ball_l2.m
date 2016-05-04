%DIST2BALL_L2 Squared distance from a Euclidean ball of given center and radius
%
%   DIST2BALL_L2(rho, c, w) builds the function
%       
%       f(x) = (w/2)*dist^2(x,B) where B is the ball ||x-c|| <= rho
%
%   If c is not provided, c = 0. If also rho is not provided, rho = 1.
%   Default weight is w = 1.
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

function obj = dist2Ball_l2(rho, c, weight)
    if nargin < 1 || isempty(rho)
        rho = 1;
    end
    if nargin < 2 || isempty(c)
        c = 0;
    end
    if nargin < 3 || isempty(weight)
        weight = 1;
    end
    if ~isscalar(weight) || weight <= 0
        error('third argument (weight) must be a positive scalar');
    end
    obj.makef = @() @(x) call_dist2Ball_l2_f(x, rho, c, weight);
    obj.L = weight;
end

function [val, grad] = call_dist2Ball_l2_f(x, rho, c, weight)
    xc = x-c;
    nxc = norm(xc);
    if nxc <= rho
        proj = x;
    else
        proj = c + (rho/nxc)*xc;
    end
    diff = x-proj;
    val = (0.5*w)*(diff'*diff);
    grad = weight*diff;
end
