function [flag, cache, cache_z, ops] = Cache_CheckGamma(cache, gam, bet)

[cache, ops] = Cache_ProxGradStep(cache, gam);
cache_z = Cache_Init(cache.prob, cache.z, gam);
[cache_z, ops1] = Cache_Evalf(cache_z);
ops = Ops_Sum(ops, ops1);
fz = cache_z.fx;
if fz <= cache.fx - cache.gradfx'*cache.FPR + (1-bet)/(2*gam)*(cache.normFPR^2) + 1e-6*abs(cache.fx)
    flag = 1;
else
    flag = 0;
end

end
