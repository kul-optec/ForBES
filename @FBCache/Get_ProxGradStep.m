function z = Get_ProxGradStep(cache, gam)

if nargin < 2
    gam = cache.gam;
end

gam0 = cache.gam;

if ~cache.flagGradStep || gam0 ~= gam
    cache.Get_GradStep(gam);
end

prob = cache.prob;

if ~cache.flagProxGradStep || gam0 ~= gam
    if prob.isthereD
      mugam = prob.mu*gam;
      [z, cache.gz] = prob.callg(prob.D*cache.y, mugam);
      cache.z = cache.y + prob.D'*(z - prob.D*cache.y)/prob.mu;
    else
      [cache.z, cache.gz] = prob.callg(cache.y, gam);
    end
    if cache.flagOps
        cache.ops.addproxg();
        cache.ops.addg();
    end
    cache.FPR = cache.x-cache.z;
    cache.normFPR = norm(cache.FPR(:));
    cache.gam = gam;
    cache.flagProxGradStep = true;
    cache.flagFBE = false;
    cache.flagGradFBE = false;
end

z = cache.z;
