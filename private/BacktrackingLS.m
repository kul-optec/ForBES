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

function [t, cachet, cachet1, ops, exitflag] = BacktrackingLS(cache, dir, t0, lsopt, ref, lin)

    if nargin < 5, ref = cache.FBE; end
    if nargin < 6, lin = 0; end

    [cache, ops] = CacheLineSearch(cache, dir);
    
    cachet1 = [];
    
    gam = cache.gam;
    
    t = t0;
    exitflag = 1;
    
    for i = 1:lsopt.nLS
        [cachet, ops1] = DirFBE(cache, t, 1);
        ops = OpsSum(ops, ops1);
        ft = cachet.FBE;
        if ft <= ref + t*lin
            exitflag = 0;
            break;
        end
        t = 0.5*t;
        if t <= lsopt.progTol
            exitflag = 2;
            break
        end
    end
    
    if exitflag == 0 && lsopt.testGamma
        [flagGamma, cachet, cachet1, ops1] = CheckGamma(cachet, gam, lsopt.beta);
        ops = OpsSum(ops, ops1);
        exitflag = flagGamma-1; % because CheckGamma returns 1 (good gamma) or 0 (bad gamma)
    end
end
