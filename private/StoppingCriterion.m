function stop = StoppingCriterion(cache, tol)
    normdiff = norm(cache.diff, 'inf');
    if normdiff <= 1e-14, stop = 1; return; end
    if normdiff/cache.gam <= tol, stop = 1; return; end
    stop = 0;
end