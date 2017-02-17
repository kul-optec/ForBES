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
%   and false otherwise.

function [flag, cache_bar] = Backtrack_Gamma(cache, bet)

flag = false;
[gamma_ok, cache_bar] = cache.Check_Gamma(bet);
while ~gamma_ok
    cache.Set_Gamma(cache.Get_Gamma()/2);
    flag = true;
    [gamma_ok, cache_bar] = cache.Check_Gamma(bet);
end
