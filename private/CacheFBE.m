function [cache, ops] = CacheFBE(cache, gam)

ops = OpsInit();

if cache.flagProxGradStep == 0 || cache.gam ~= gam
    [cache, ops] = CacheProxGradStep(cache, gam);
end

cache.diff = cache.z-cache.x;
cache.normdiff = norm(cache.diff);
cache.FBE = cache.fx + cache.gz + cache.gradfx'*cache.diff + (0.5/gam)*(cache.normdiff^2);
cache.gam = gam;

cache.flagFBE = 1;
