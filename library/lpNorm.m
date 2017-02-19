%LPNORM Allocates the L1 norm function.
%
%   LPNORM(type, mu) builds the function
%
%       g(x) = mu*sum(|x_i|^p)
% 
%   where if type = 'onehalf' then p = 1/2,
%                   'twothirds' then p = 2/3.
%   If mu is not provided, then mu = 1.0.

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

function obj = lpNorm(type, mu)
    if nargin < 1
        error('must provide type = ''onehalf'' or ''twothirds''');
    end
    if nargin < 2
        mu = 1;
    end
    obj.isConvex = 0;
    obj.isQuadratic = 0;
    obj.isConjQuadratic = 0;
    switch type
        case 'onehalf'
            obj.makeprox = @() @(x, gam) call_l_1_2_prox(x, gam, mu);
        case 'twothirds'
            error('not yet implemented');
%             obj.makeprox = @() @(x, gam) call_l_2_3_prox(x, gam, mu);
    end
end

% The implementation of the proximal mappings is based on:
% Cao, Sun, Xu, "Fast image deconvolution using closed-form
% thresholding formulas of L_q (q=1/2,2/3) regularization" (2013)

function [prox, g] = call_l_1_2_prox(x, gam, mu)
    lam = 2*gam*mu;
    p = nthroot(54, 3)/4*nthroot(lam, 3)^2;
    absx = abs(x);
    phi = acos((lam/8)*sqrt((absx/3).^(-3)));
    ind0 = (absx <= p);
    prox(ind0, 1) = 0;
    prox(~ind0, 1) = (2/3)*(sign(x(~ind0)).*absx(~ind0).*(1+cos((2/3)*(pi-phi(~ind0)))));
    g = mu*sum(abs(prox).^(0.5));
end

% The following doesn't seem to work; must check what's wrong

function [prox, g] = call_l_2_3_prox(x, gam, mu)
    lam = 2*gam*mu;
    p = (2/3)*(3*lam^3)^(0.25);
    absx = abs(x);
    phi = sqrt(lam^(-3))*(27/16)*x.^2;
    A = (2/sqrt(3))*lam^(0.25)*sqrt(cosh(phi/3));
    ind0 = (absx <= p);
    prox(ind0, 1) = 0;
    prox(~ind0, 1) = sign(x(~ind0)).*((A(~ind0) + sqrt(2*absx(~ind0)./A(~ind0) - A(~ind0).^2))/2).^3;
    g = mu*sum(abs(prox).^(2/3));
end
