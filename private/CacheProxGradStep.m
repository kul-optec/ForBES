function [cache, ops] = CacheProxGradStep(cache, gam)

if nargin < 2
    gam = cache.gam;
end

gam0 = cache.gam;

if cache.flagGradStep == 0 || gam0 ~= gam
    [cache, ops] = CacheGradStep(cache, gam);
else
    ops = OpsInit();
end

if cache.flagProxGradStep == 0 || gam0 ~= gam
    [cache.z, cache.gz] = cache.prob.callg(cache.y, gam);
    ops.proxg = ops.proxg + 1;
    ops.g = ops.g + 1;
end

cache.diff = cache.z-cache.x;
cache.normdiff = norm(cache.diff);

cache.gam = gam;

cache.flagProxGradStep = 1;
