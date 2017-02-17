% GET_GRADSTEP(cache)
%
%   Returns the gradient step y = x - gam*gradf(x), where x and gam are
%   stored in cache.

function y = Get_GradStep(cache)

if cache.flagGradStep == true
    y = cache.y;
    return;
end

if cache.flagGradf == false
    cache.Get_Gradf();
end

gam = cache.gam;

cache.y = cache.x - gam*cache.Get_Gradf();

cache.flagGradStep = true;
y = cache.y;
