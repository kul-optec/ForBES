function stop = StoppingCriterion(cache, tol, cache_prev)

% Absolute stopping criterion on the fixed point residual
normInfFPR = norm(cache.FPR, 'inf')/cache.gam;
if normInfFPR <= tol
    stop = true;
else
    stop = false;
end

% % From sec. 8.2.3.2 of Gill, Murray, Wright (1982).
% if nargin < 3, cache_prev = []; end
% normInfFPR = norm(cache.FPR, 'inf')/cache.gam;
% absFBE = abs(cache.FBE);
% stop = ( normInfFPR <= 10*sqrt(eps) ) || ...
%        ( ~isempty(cache_prev) && normInfFPR <= nthroot(tol, 3)*(1+absFBE) && ...
%          norm(cache_prev.x-cache.x, inf) < sqrt(tol)*(1+norm(cache.x, inf)) && ...
%          cache_prev.FBE-cache.FBE < tol*(1+absFBE) )
