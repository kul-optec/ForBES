function stop = StoppingCriterion(cache, tol)
    normInfFPR = norm(cache.FPR, 'inf');
    if normInfFPR <= 1e-14, stop = 1; return; end
    if normInfFPR/cache.gam <= tol, stop = 1; return; end
    stop = 0;
end