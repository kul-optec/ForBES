function FPR = Get_FPR(cache)

if cache.flagProxGradStep == true
    FPR = cache.FPR;
    return;
end

cache.Get_ProxGradStep();
FPR = cache.FPR;
