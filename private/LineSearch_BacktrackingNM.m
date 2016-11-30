function [t, cachet, cachet1, ops, lsopt, exitflag] = LineSearch_BacktrackingNM(cache, dir, slope, t0, lsopt, it, restart, ref, lin, const)

if it == 1 || restart
    lsopt.Q = 1;
    lsopt.C = ref;
else
    newQ = lsopt.eta*lsopt.Q+1;
    lsopt.C = (lsopt.eta*lsopt.Q*lsopt.C + ref)/newQ;
    lsopt.Q = newQ;
end

[t, cachet, cachet1, ops, lsopt, exitflag] = LineSearch_Backtracking(cache, dir, slope, t0, lsopt, it, restart, lsopt.C, lin, const);

end
