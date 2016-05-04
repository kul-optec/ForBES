%HINGELOSS Allocates the hinge loss function.
%
%   HINGELOSS(mu, b) builds the function
%       
%       g(x) = mu*sum(max(0, 1-b.*x))
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

function obj = hingeLoss(mu, b)
    %
    % Only the proximal mapping is available for this function
    %
    obj.makeprox = @() @(x, gam) call_hingeLoss_prox(x, gam, mu, b);
end

function [prox, g] = call_hingeLoss_prox(x, gam, mu, b)
    bx = b.*x; ind = bx < 1;
    prox(ind,1) = b(ind).*min(bx(ind)+gam*mu,1);
    prox(~ind,1) = x(~ind);
    g = mu*sum(max(0,1-b.*prox));
end
