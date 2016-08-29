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

function cache = CacheInit(prob, x, gam)

cache.prob = prob;
cache.x = x;
cache.gam = gam;

cache.flagInit = 1;
cache.flagEvalf = 0;
cache.flagGradStep = 0;
cache.flagProxGradStep = 0;
cache.flagFBE = 0;
cache.flagGradFBE = 0;
cache.flagLineSearch = 0;
