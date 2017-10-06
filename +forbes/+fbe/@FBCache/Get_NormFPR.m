% GET_NORMFPR(cache)
%
%   Returns the Euclidean norm of the fixed-point residual vector x-z,
%   where x is the point to which cache refers, and z is the proximal gradient
%   point.

function normFPR = Get_NormFPR(cache)

if cache.flagProxGradStep == true
    normFPR = cache.normFPR;
    return;
end

cache.Get_ProxGradStep();
normFPR = cache.normFPR;
