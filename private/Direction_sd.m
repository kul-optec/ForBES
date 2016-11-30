% Steepest descent

function [dir, tau0, cache] = Direction_sd(prob, opt, it, restart, sk, yk, v, cache)

v = v(:);

dir = -v;
tau0 = 1.0/norm(dir, inf);

end
