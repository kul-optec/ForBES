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

function [cachet, t, cnt, exitflag] = ArmijoLS(prob, gam, cache, df0, t0, lsopt, fr)
%ARMIJOLS - computes a steplength t > 0 so that it satisfies the Armijo condition
%
% f(t) <= f(0) + delta*f'(0)
%
% exitflag = -1: gam is not small enough
% exitflag =  0: acceptable steplength was found
% exitflag =  1: maximum number of backtracking steps exceeded
% exitflag =  2: no further progress can be made

    %      Q, A, C, f2, proxg
    cnt = [0, 0, 0, 0, 0];
    arm_hi = lsopt.delta*df0;
    t = t0;
    exitflag = 1;
    if nargin >= 7
        f0 = fr;
    else
        f0 = cache.FBE;
    end
    for i = 1:lsopt.nLS
        [cachet, cnt1] = DirFBE(prob, gam, t, cache, 1);
        cnt = cnt+cnt1;
        ft = cachet.FBE;
        if ft <= f0 + t*arm_hi
            exitflag = 0;
            break;
        end
        if i == 1 %quadratic interpolation
            tn = ArmijoQuadInterp(f0,df0,t,ft);
        else %cubic interpolation
            tn = ArmijoCubInterp(f0,df0,told,ftold,t,ft);
        end
        if tn <= 0
            tn = 0.5*t;
        end
        told = t;ftold = ft;
        t = tn;
        if t <= lsopt.progTol
            exitflag = 2;
            break
        end
    end
    if exitflag == 0 && lsopt.testGamma
        [fz, cnt1] = Evaluatef(prob, cachet.z);
        cnt = cnt+cnt1;
        % check whether gam is small enough
        if fz + cachet.gz - cachet.FBE > 1e-14*(1+abs(cachet.FBE))
            exitflag = -1;
        end
    end
end

function t = ArmijoQuadInterp(f0,df0,t,ft)
    % Minimizer of interpolant belongs to [0,t1]
    tdf0 = t*df0;
    q = ft-f0-tdf0;
    if q > 0%quadratic is strongly convex
        t = -(tdf0*t)/(2*q);
    else
        t = -1;
    end
end

function t = ArmijoCubInterp(f,df,t0,f0,t1,f1)
    % Minimizer of interpolant belongs to [0,t1]
    t02 = t0^2;
    t12 = t1^2;
    ab = 1/(t02*t12*(t1-t0))*[t02 -t12;-t0^3 t1^3]*[f1-f-df*t1;f0-f-df*t0];
    a = ab(1);
    b = ab(2);
    t = (-b+sqrt(b^2-3*a*df))/(3*a);
end
