function [t, cachet, cachet1, ops, lsopt, exitflag] = LineSearch_Backtracking...
    (cache, direction, slope, t0, lsopt, it, restart, ref, lin, const)

    if nargin < 8, ref = cache.FBE; end
    if nargin < 9, lin = 0.0; end
    if nargin < 10, const = 0.0; end

    [cache, ops] = Cache_LineSearch(cache, direction);

    cachet1 = [];

    gam = cache.gam;

    t = t0;
    exitflag = 1;

    for i = 1:lsopt.nLS
        [cachet, ops1] = Cache_LineFBE(cache, t, 1);
        ops = Ops_Sum(ops, ops1);
        ft = cachet.FBE;
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

    if exitflag == 0 && lsopt.testGamma
        [isGammaOK, cachet, cachet1, ops1] = Cache_CheckGamma(cachet, gam, lsopt.beta);
        ops = Ops_Sum(ops, ops1);
        exitflag = isGammaOK-1; % because CheckGamma returns 1 (good gamma) or 0 (bad gamma)
    end
end
