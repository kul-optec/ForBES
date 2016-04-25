%MOREAUENVELOPE Defines the Moreau envelope of a proximable function
%
%   MOREAUENVELOPE(f, gam) where f is a function object, and gam > 0 is a
%   real number, returns function F(x) defined as
%       
%       F(x) = f(prox(x, gam)) + 1/(2*gam)*||x - prox(x, gam)||^2
%
%   where prox is the proximal mapping associated with f.
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

function obj = moreauEnvelope(obj1, gam)
    if ~isfield(obj1, 'makeprox')
        error('the function does not have the proximal mapping defined');
    end
    if gam <= 0
        error('second argument must be a positive real number');
    end
    proxf1 = obj1.makeprox();
    obj.L = 1/gam;
    obj.makef = @() @(x) call_moreauEnvelope_f1(x, proxf1, gam);
    if isfield(obj1, 'isConvex')
        obj.isConvex = obj1.isConvex;
    end
end

function [v, grad] = call_moreauEnvelope_f1(x, proxf1, gam)
    [z, f1z] = proxf1(x, gam);
    res = x - z;
    v = f1z + 0.5/gam*norm(res, 2)^2;
    if nargout >= 2
        grad = res/gam;
    end
end
