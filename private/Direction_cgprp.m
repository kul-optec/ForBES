% CG method of Polak-Ribiere-Polyak

function [dir, tau0, cache] = Direction_cgprp(prob, opt, it, restart, sk, yk, v, cache)

sk = sk(:);
yk = yk(:);
v = v(:);

if it == 1 || restart
    dir = -v; % Initially use steepest descent direction
    tau0 = 1.0/norm(dir, inf);
else
    beta = max((v'*yk)/cache.sqnorm_prev_v,0);
    dir = -v + beta*cache.prev_dir;
    if dir'*v >= 0 % restart if not descent direction
        dir = -v;
        tau0 = 1.0/norm(dir, inf);
        cache.cntSkip = cache.cntSkip+1;
    else
        tau0 = (sk'*sk)/(sk'*yk);
    end
end

cache.sqnorm_prev_v = norm(v)^2;
cache.prev_dir = dir;

end
