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

function [cache, ops] = CacheLineSearch(cache, dir)

if cache.flagLineSearch == 1
    if norm(dir-cache.dir) == 0
        ops = OpsInit();
        return;
    end
end

if cache.flagGradStep == 0
    [cache, ops] = CacheGradStep(cache, cache.gam);
else
    ops = OpsInit();
end

prob = cache.prob;
cache.dir = dir;

if prob.istheref1
    if prob.isthereC1
        cache.C1dir = prob.C1*dir;
        cache.QC1dir = prob.Q(cache.C1dir);
        cache.C1tQC1dir = prob.C1'*cache.QC1dir;
        ops.C1 = ops.C1 + 2;
    else
        cache.C1dir = dir;
        cache.QC1dir = prob.Q(cache.C1dir);
        cache.C1tQC1dir = cache.QC1dir;
    end
    ops.gradf1 = ops.gradf1 + 1;
    cache.f1linear = cache.gradf1x(:)'*dir(:);
    cache.f1quad = cache.C1dir(:)'*cache.QC1dir(:);
end

if prob.istheref2
    if prob.isthereC2
        cache.C2dir = prob.C2*dir;
        ops.C2 = ops.C2 + 1;
    else
        cache.C2dir = dir;
    end
end

if prob.istherelin
    cache.lindir = prob.l(:)'*dir(:);
end

cache.flagLineSearch = 1;
