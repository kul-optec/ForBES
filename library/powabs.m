%POWABS POWer of ABSolute value function.
%
%   POWABS(p, c) builds the function
%       
%       g(x) = c*|x|^p
%   
%   If c is not provided, it is assumed to be 1. If also p is not provided,
%   it is assumed to be 1.
%
%   See Combettes, Dung, Vu, Dualization of Signal Recovery Problems, Example 2.15
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

function obj = powabs(p, c)
    if nargin<2 || isempty(c)
        c = 1;
        if nargin<1 || isempty(p)
            p = 1;
        end
    end
    pows = [1;4/3;3/2;2;3;4];
    if all((p == pows) == 0)
        error('Prox is not computable')
    end
    obj.makeprox = @() @(x, gam) call_powabs_prox(x, gam, p, c);
end

function [prox, val] = call_powabs_prox(x, gam, p, c)
    gam = c*gam;
    switch p
        case 1
            prox = sign(x)*max(abs(x)-gam,0);
        case 4/3
            rho = sqrt(x^2+256*gam^3/729);
            prox = x + (4*gam)/(3*2^(1/3))*(abs(rho-x)^(1/3)-abs(rho+x)^(1/3));
        case 3/2
            prox = x + 9*gam^2*sign(x)*(1-sqrt(1+16*abs(x)/(9*gam^2)))/8;
        case  2
            prox = x/(1+2*gam);
        case  3
            prox = sign(x)*(sqrt(1+12*gam*abs(x))-1)/(6*gam);
        case  4
            rho = sqrt(x^2+1/(27*gam));
            prox = abs((rho + x)/(8*gam))^(1/3)-abs((rho - x)/(8*gam))^(1/3);
    end
    val = c*abs(prox)^p;
end
