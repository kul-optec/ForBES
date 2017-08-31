% GET_PROXGRADSTEP(cache)
%
%   Returns the proximal gradient point z = prox_{gam*g}(x - gam*gradf(x)),
%   where x and gam are stored in cache.

function [z, gz] = Get_ProxGradStep(cache)

if cache.flagProxGradStep == true
    z = cache.z;
    gz = cache.gz;
    return;
end

if cache.flagGradStep == false
    cache.Get_GradStep();
end

prob = cache.prob;
gam = cache.gam;

if prob.isthereD
    mugam = prob.mu*gam;
    [z, cache.gz] = prob.callg(prob.D*cache.y, mugam);
    cache.z = cache.y + prob.D'*(z - prob.D*cache.y)/prob.mu;
else
    [cache.z, cache.gz] = prob.g.prox(cache.y, gam);
end
if cache.flagOps
    cache.ops.addproxg();
    cache.ops.addg();
end
cache.FPR = cache.x-cache.z;
cache.normFPR = norm(cache.FPR(:));

cache.flagProxGradStep = true;
z = cache.z;
gz = cache.gz;
