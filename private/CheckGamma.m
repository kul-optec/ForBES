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

function [flag, cache, cache_z, ops] = CheckGamma(cache, gam, bet)

[cache, ops] = CacheProxGradStep(cache, gam);
cache_z = CacheInit(cache.prob, cache.z, gam);
[cache_z, ops1] = CacheEvalf(cache_z);
ops = OpsSum(ops, ops1);
fz = cache_z.fx;
if fz <= cache.fx - cache.gradfx'*cache.FPR + (1-bet)/(2*gam)*(cache.normFPR^2) + 1e-6*abs(cache.fx)
    flag = 1;
else
    flag = 0;
end
