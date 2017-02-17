% BACKTRACK_GAMMA(cache, bet)
%
%   Reduces the stepsize gam in cache until the following inequality holds
%
%       f(z) <= f(x) - gradf(x)'*(x - z) + (1-bet)/(2*gam)||x - z||^2
%
%   where x is the point to which cache refers and z is the proximal gradient
%   step with stepsize gam.
%
%   Returns true if gamma has changed (i.e., some backtracks have occurred)
%   and false otherwise. The second return argument cache_pg is the cache at
%   the proximal gradient point z (which contains some information, and
%   therefore may be re-used).

function [flag, cache_pg] = Backtrack_Gamma(cache, bet)

flag = false;
[gamma_ok, cache_pg] = cache.Check_Gamma(bet);
while ~gamma_ok
    cache.Set_Gamma(cache.Get_Gamma()/2);
    flag = true;
    [gamma_ok, cache_pg] = cache.Check_Gamma(bet);
end
