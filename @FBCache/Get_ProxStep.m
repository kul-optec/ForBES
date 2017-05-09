% GET_PROXSTEP(cache, y)
%
%   Returns the proximal point z = prox_{gam*g}(y), where y is some given point.

function [z, gz] = Get_ProxStep(cache, y)

prob = cache.prob;
gam = cache.gam;

if prob.isthereD
    mugam = prob.mu*gam;
    [z, gz] = prob.callg(prob.D*y, mugam);
    z = y + prob.D'*(z - prob.D*y)/prob.mu;
else
    [z, gz] = prob.callg(y, gam);
end
if cache.flagOps
    cache.ops.addproxg();
    cache.ops.addg();
end
