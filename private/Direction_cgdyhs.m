% CG, hybrid method Dai-Yuan/Hestenes-Stiefel

function [dir, tau0, cache] = Direction_cgdyhs(prob, opt, it, restart, sk, yk, v, cache)

sk = sk(:);
yk = yk(:);
v = v(:);

if it == 1 || restart
    dir = -v; % Initially use steepest descent direction
    tau0 = 1.0/norm(dir, inf);
else
    betaDY = (v'*v)/(sk'*yk);
    betaHS = (v'*yk)/(sk'*yk);
    beta = max(0.0,min(betaHS,betaDY));
    dir = -v + beta*sk;
    if dir'*v >= 0 % restart if not descent direction
        dir = -v;
        tau0 = 1.0/norm(dir, inf);
        cache.cntSkip = cache.cntSkip+1;
    else
        tau0 = (sk'*sk)/(sk'*yk);
    end
end

end
