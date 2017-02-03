function FPR = Get_FPR(cache, gam)

if nargin < 2
    gam = cache.gam;
end

gam0 = cache.gam;

if ~cache.flagProxGradStep || gam0 ~= gam
    cache.Get_ProxGradStep(gam);
end

FPR = cache.FPR;
