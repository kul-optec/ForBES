% GET_FBE(cache)
%
%   Returns the value of the FBE (the Forward-Backward envelope) at the point x
%   to which cache refers.

function FBE = Get_FBE(cache)

if cache.flagFBE == true
    FBE = cache.FBE;
    return;
end

if cache.flagProxGradStep == false
    cache.Get_ProxGradStep();
end

gam = cache.gam;

cache.FBE = cache.fx + cache.gz - cache.gradfx(:)'*cache.FPR(:) + (0.5/gam)*(cache.normFPR^2);

cache.flagFBE = true;
FBE = cache.FBE;
