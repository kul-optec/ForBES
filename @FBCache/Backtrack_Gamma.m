function [flag, cache_bar] = Backtrack_Gamma(cache, bet)

flag = false;
[gamma_ok, cache_bar] = cache.Check_Gamma(bet);
while ~gamma_ok
    cache.gam = cache.gam/2;
    flag = true;
    [gamma_ok, cache_bar] = cache.Check_Gamma(bet);
end
