%L2NORM Allocates the L2 norm function.
%
%   L2NORM(mu) builds the function
%       
%       g(x) = mu*||x||_2
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

function obj = l2Norm(mu)
    if nargin < 1
        mu = 1;
    end
    obj.makeprox = @() @(x, gam) call_l2Norm_prox(x, gam, mu);
end

function [prox, val] = call_l2Norm_prox(x, gam, mu)
    normx = sqrt(x'*x);
    mugam = mu*gam;
    if normx <= mugam
        prox = zeros(length(x),1);
        val = 0;
    else
        scal = (1-mugam/normx);
        prox = (1-mugam/normx)*x;
        val = mu*scal*normx;
    end
end
