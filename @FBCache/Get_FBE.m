function FBE = Get_FBE(cache, gam)

if nargin < 2
    gam = cache.gam;
end

gam0 = cache.gam;

if ~cache.flagProxGradStep || gam0 ~= gam
    cache.Get_ProxGradStep(gam);
end

if ~cache.flagFBE || gam0 ~= gam
    cache.FBE = cache.fx + cache.gz - cache.gradfx(:)'*cache.FPR(:) + (0.5/gam)*(cache.normFPR^2);
    cache.gam = gam;
    cache.flagFBE = true;
end

FBE = cache.FBE;
