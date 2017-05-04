% GET_G(cache)
%
%   Returns the value of g, the nonsmooth term of the problem, evaluated at
%   z, the proximal gradient point to which cache refers.

function gz = Get_g(cache)

if cache.flagProxGradStep
    gz = cache.gz;
    return;
end

cache.Get_ProxGradStep();
gz = cache.gz;
