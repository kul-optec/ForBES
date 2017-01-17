function [t, cachet, cachet1, lsopt, exitflag] = LineSearch_BacktrackingArmijo(cache, dir, slope, t0, lsopt, adaptive, it, restart, varargin)

ref = cache.FBE; % f(0)
lin = lsopt.delta*slope; % delta f'(0)
[t, cachet, cachet1, lsopt, exitflag] = LineSearch_Backtracking(cache, dir, slope, t0, lsopt, adaptive, it, restart, ref, lin);

end
