function [t, cachet, cachet1, lsopt, exitflag] = LineSearch_Backtracking...
    (cache, direction, slope, t0, lsopt, adaptive, it, restart, ref, lin, const)

    if nargin < 9, ref = cache.Get_FBE(); end
    if nargin < 10, lin = 0.0; end
    if nargin < 11, const = 0.0; end

    cache.Set_Directions(direction);

    cachet1 = [];

    t = t0;
    exitflag = 1;

    for i = 1:lsopt.nLS
        cachet = cache.Get_CacheLine(t, 1);
        ft = cachet.Get_FBE();
        if ft <= ref + t*lin + 1e-14*abs(ref) + const
            exitflag = 0;
            break;
        end
        t = 0.5*t;
        if t <= lsopt.progTol
            exitflag = 2;
            break
        end
    end

    if exitflag == 0 && adaptive
        [flag, cachet1] = cachet.Backtrack_Gamma(lsopt.beta);
        exitflag = -1*flag;
    end
end
