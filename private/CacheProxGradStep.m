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

function [cache, ops] = CacheProxGradStep(cache, gam)

if nargin < 2
    gam = cache.gam;
end

gam0 = cache.gam;

if cache.flagGradStep == 0 || gam0 ~= gam
    [cache, ops] = CacheGradStep(cache, gam);
else
    ops = OpsInit();
end

if cache.flagProxGradStep == 0 || gam0 ~= gam
    [cache.z, cache.gz] = cache.prob.callg(cache.y, gam);
    ops.proxg = ops.proxg + 1;
    ops.g = ops.g + 1;
    cache.FPR = cache.x-cache.z;
    cache.normFPR = norm(cache.FPR(:));
    cache.gam = gam;
    cache.flagProxGradStep = 1;
    cache.flagFBE = 0;
    cache.flagGradFBE = 0;
end
