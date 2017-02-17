% GET_FPR(cache)
%
%   Returns the fixed-point residual vector x-z, where x is the point to which
%   cache refers, and z is the proximal gradient point.

function FPR = Get_FPR(cache)

if cache.flagProxGradStep == true
    FPR = cache.FPR;
    return;
end

cache.Get_ProxGradStep();
FPR = cache.FPR;
