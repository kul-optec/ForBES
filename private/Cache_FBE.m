function [cache, ops] = Cache_FBE(cache, gam)

ops = Ops_Init();

gam0 = cache.gam;

if cache.flagProxGradStep == 0 || gam0 ~= gam
    [cache, ops] = Cache_ProxGradStep(cache, gam);
end

if cache.flagFBE == 0 || gam0 ~= gam
    cache.FBE = cache.fx + cache.gz - cache.gradfx(:)'*cache.FPR(:) + (0.5/gam)*(cache.normFPR^2);
    cache.gam = gam;
    cache.flagFBE = 1;
end

end
