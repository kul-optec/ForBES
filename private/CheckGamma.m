function [flag, cache, cache_z, ops] = CheckGamma(cache, gam, bet)

[cache, ops] = CacheProxGradStep(cache, gam);
cache_z = CacheInit(cache.prob, cache.z, gam);
[cache_z, ops1] = CacheEvalf(cache_z);
ops = OpsSum(ops, ops1);
fz = cache_z.fx;
if fz <= cache.fx + cache.gradfx'*cache.diff + 0.5*gam*(1-bet)*(cache.normdiff^2)
    flag = 1;
else
    flag = 0;
end
