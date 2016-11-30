function [t, cachet, cachet1, ops, lsopt, exitflag] = LineSearch_BacktrackingArmijo(cache, dir, slope, t0, lsopt, it, restart, varargin)

ref = cache.FBE; % f(0)
lin = lsopt.delta*slope; % delta f'(0)
[t, cachet, cachet1, ops, lsopt, exitflag] = LineSearch_Backtracking(cache, dir, slope, t0, lsopt, it, restart, ref, lin);

end
