function normFPR = Get_NormFPR(cache)

if cache.flagProxGradStep == true
    normFPR = cache.normFPR;
    return;
end

cache.Get_ProxGradStep();
normFPR = cache.normFPR;
