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

function obj = distBall_l2(rho,c,weight)
    % Proximal mapping of w*dist(x,C) where C is the l2 ball ||x-c||<=rho
    if nargin<3 || isempty(weight)
        weight = 1;
    end
    if nargin<2 || isempty(c)
        c = 0;
    end
    if nargin<1 || isempty(rho)
        rho = 1;
    end

    obj.makeprox = @() @(x, gam) call_distBall_l2_prox(x, gam, rho, c, weight);
end

function [prox, val] = call_distBall_l2_prox(x, gam, rho, c, weight)
    gam = weight*gam;
    xc = x-c;
    nxc = norm(xc);
    if nxc <= rho
        prox = x;
        val = 0;
    else
        scale = rho/nxc;
        if nxc > gam/(1-scale);
            prox = x-(gam/nxc)*xc;
            xc = prox-c;
            nxc = norm(xc);
            if nxc <= rho
                val = 0;
            else
                val = weight*(1-rho/nxc)*nxc;
            end
        else
            prox = c + scale*xc;
            val = 0;
        end
    end
end
