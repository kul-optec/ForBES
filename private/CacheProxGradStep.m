function [cache, ops] = CacheProxGradStep(cache, gam)

ops = OpsInit();
gam0 = cache.gam;

if cache.flagGradStep == 0 || gam0 ~= gam
    [cache, ops] = CacheGradStep(cache, gam);
end

if cache.flagProxGradStep == 0 || gam0 ~= gam
    [cache.z, cache.gz] = cache.prob.callg(cache.y, gam);
    ops.proxg = ops.proxg + 1;
    ops.g = ops.g + 1;
end

cache.gam = gam;

cache.flagProxGradStep = 1;
