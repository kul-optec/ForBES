% INDEPISQUAREDNORM Indicator of the epigraph of the squared norm
%
%  INDEPISQUAREDNORM builds the indicator of the epigraph of the squared
%  norm, that is, the set
%
%   C = {(x,t) : ||x||^2 <= t}
%


% Copyright (C) 2015-2017, KUL-FobBES (https://github.com/kul-forbes)
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

function obj = indEpiSquaredNorm()
obj.isConvex = 1;
obj.makeprox = @() @(x, gam) epipr_sqnorm(x(1:end-1), x(end));

function [prox, val] = epipr_sqnorm(x,z)
if (x'*x <= z)
    prox =[x; z]; val = 0;
    return;
end
theta = 1 - 2 * z;
r = cubic_roots(theta, x);
for i=1:length(r),   % Pick the right root
    x_ = x/(1 + 2*(r(i) - z));
    if abs(norm(x_)^2-r(i)) < 1e-6,
        z_ = r(i);
        break;
    end
end
prox=[x_; z_];
val = 0;

function [r, status] = cubic_roots(theta, x)
b=4*theta; c=theta^2; d=-x'*x;
D  = 72*b*c*d - 4*b^3*d +b^2*c^2 - 16*c^3 - 432*d^2;
D0 = b^2 - 12*c;
status.D = D; status.D0 = D0;
if abs(D)<1e-14,
    if abs(D0)<1e-14, % one triple root
        r = -b/12;
        status.msg = 'one triple root';
    else % a double root and a single one
        r = zeros(2,1);
        r(1) = (16*b*c - 144*d - b^3)/(4*D0); % single
        r(2) = (36*d - b*c)/(2*D0); % double (cannot be)
        status.msg = 'double plus single';
    end
    return;
end
r = roots([4 b c d]); % eigenvalues of matrix