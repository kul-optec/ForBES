%HINGELOSS Allocates the hinge loss function.
%
%   HINGELOSS(b, mu) builds the function
%       
%       g(x) = sum(mu.*max(0, 1-b.*x))
%
%   If mu is not provided then mu = 1.
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

function obj = hingeLoss(b, mu)
    if nargin < 1 || isempty(b), error('first argument b is required'); end
    if nargin < 2 || isempty(mu), mu = 1; end
    if any(mu < 0)
        error('second argument mu must be nonnegative');
    end
    obj.makeprox = @() @(x, gam) call_hingeLoss_prox(x, gam, mu, b);
    obj.isConvex = 1;
end

function [prox, g] = call_hingeLoss_prox(x, gam, mu, b)
    bx = b.*x; ind = bx < 1;
    prox(ind,1) = b(ind).*min(bx(ind)+gam*mu,1);
    prox(~ind,1) = x(~ind);
    g = sum(mu.*max(0,1-b.*prox));
end
