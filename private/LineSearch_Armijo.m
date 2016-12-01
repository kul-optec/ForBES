function [t, cachet, cachet1, ops, exitflag] = LineSearch_Armijo(cache, dir, slope, t0, lsopt, ref)
%ARMIJOLS - computes a steplength t > 0 so that it satisfies the Armijo condition
%
% f(t) <= f(0) + delta*f'(0)
%
% exitflag = -1: gam is not small enough
% exitflag =  0: acceptable steplength was found
% exitflag =  1: maximum number of backtracking steps exceeded
% exitflag =  2: no further progress can be made

    % precompute stuff for the line search
    [cache, ops] = CacheLineSearch(cache, dir);

    cachet1 = [];

    gam = cache.gam;

    arm_hi = lsopt.delta*slope;
    t = t0;
    exitflag = 1;
    if nargin >= 6
        f0 = ref;
    else
        f0 = cache.FBE;
    end
    for i = 1:lsopt.nLS
        [cachet, ops1] = LineFBE(cache, t, 1);
        ops = OpsSum(ops, ops1);
        ft = cachet.FBE;
        if ft <= f0 + t*arm_hi
            exitflag = 0;
            break;
        end
        if i == 1 %quadratic interpolation
            tn = ArmijoQuadInterp(f0, slope, t, ft);
        else %cubic interpolation
            tn = ArmijoCubInterp(f0, slope, told, ftold, t, ft);
        end
        if tn <= 0
            tn = 0.5*t;
        end
        told = t;
        ftold = ft;
        t = tn;
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
