function Set_Gamma(cache, gam)

cache.gam = gam;

% When gamma changes, the gradient step, proximal-gradient step
% and value/gradient of the FBE need to be recomputed.
% Reset flags for this purpose.

cache.flagGradStep = false;
cache.flagProxGradStep = false;
cache.flagFBE = false;
cache.flagGradFBE = false;
