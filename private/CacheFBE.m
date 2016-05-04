function [cache, ops] = CacheFBE(cache, gam)

ops = OpsInit();

gam0 = cache.gam;

if cache.flagProxGradStep == 0 || gam0 ~= gam
    [cache, ops] = CacheProxGradStep(cache, gam);
end

if cache.flagFBE == 0 || gam0 ~= gam
    cache.FBE = cache.fx + cache.gz + cache.gradfx'*cache.diff + (0.5/gam)*(cache.normdiff^2);
    cache.gam = gam;
    cache.flagFBE = 1;
end
