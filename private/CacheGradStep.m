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

function [cache, ops] = CacheGradStep(cache, gam)

ops = OpsInit();

if nargin < 2
    gam = cache.gam;
end

if cache.flagGradStep == 1
    if cache.gam ~= gam
        cache.gam = gam;
        cache.y = cache.x - gam*cache.gradfx;
        cache.flagProxGradStep = 0;
        cache.flagFBE = 0;
        cache.flagGradFBE = 0;
    end
    return;
end

if cache.flagEvalf == 0
    [cache, ops] = CacheEvalf(cache);
end

cache.gam = gam;
prob = cache.prob;

if prob.istheref1
    if prob.isthereC1
        cache.gradf1x = prob.C1'*cache.gradf1res1x;
        ops.C1 = ops.C1 + 1;
    end
else
    cache.gradf1x = 0;
end

if prob.istheref2
    if prob.isthereC2
        if prob.useHessian
            [~, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
        else
            [~, gradf2res2x] = prob.callf2(cache.res2x);
            cache.gradf2res2x = gradf2res2x;
        end
        cache.gradf2x = prob.C2'*gradf2res2x;
        ops.C2 = ops.C2 + 1;
    else
        if prob.useHessian
            [~, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
        else
            [~, gradf2res2x] = prob.callf2(cache.res2x);
            cache.gradf2res2x = gradf2res2x;
        end
        cache.gradf2x = gradf2res2x;
    end
    ops.gradf2 = ops.gradf2 + 1;
else
    cache.gradf2x = 0;
end

if prob.istherelin
    cache.gradfx = cache.gradf1x + cache.gradf2x + prob.l;
else
    cache.gradfx = cache.gradf1x + cache.gradf2x;
end

cache.y = cache.x - gam*cache.gradfx;

cache.flagGradStep = 1;
