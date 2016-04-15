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

function [t, cachet, cachet1, ops, exitflag] = BacktrackingLS(cache, t0, lsopt)

    ops = OpsInit();
    cachet1 = [];
    
    prob = cache.prob;
    gam = cache.gam;
    
    t = t0;
    exitflag = 1;
    f0 = cache.FBE;
    
    for i = 1:lsopt.nLS
        [cachet, ops1] = DirFBE(cache, t, 1);
        ops = OpsSum(ops, ops1);
        ft = cachet.FBE;
        if ft <= f0
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
        cachet1 = CacheInit(prob, cachet.z, gam);
        [cachet1, ops1] = CacheEvalf(cachet1);
        ops = OpsSum(ops, ops1);
        fz = cachet1.fx;
        % check whether gam is small enough
        if fz + cachet.gz > cachet.FBE %+ 1e-14*(1+abs(cachet.FBE))
            exitflag = -1;
        end
    end
end
