% CHECK_GAMMA(cache, bet)
%
%   Checks whether the following inequality is satisfied by cache
%
%       f(z) <= f(x) - gradf(x)'*(x - z) + (1-bet)/(2*gam)||x - z||^2
%
%   where x is the point to which cache refers, z is the proximal gradient
%   step with stepsize gam, and gam is the stepsize parameter in cache.
%
%   Returns true if the inequality holds, false otherwise. The second return
%   argument cache_pg is the cache at the proximal gradient point z (which
%   contains some information, and therefore may be re-used).
%
%   If the inequality holds, then gam <= (1-bet)/L, where L is a "local
%   Lipschitz constant" of gradf, i.e., such that the quadratic upper bound
%   for f holds between the points x and z.

function [flag, cache_pg] = Check_Gamma(cache, bet)

z = cache.Get_ProxGradStep();
cache_pg = FBCache(cache.prob, z, cache.gam, cache.ops);
fz = cache_pg.Get_f();

% We should do the following, but maybe accessing the fields directly has
% smaller overhead (no function call, no flag checks, and so on)

% fx = cache.Get_f();
% gradfx = cache.Get_Gradf();
% fprx = cache.Get_FPR();

% So we directly do the following instead

if fz <= cache.fx - cache.gradfx(:)'*cache.FPR(:) + (1-bet)/(2*cache.gam)*(cache.normFPR()^2) + 1e-6*abs(cache.fx)
    flag = true;
else
    flag = false;
end
