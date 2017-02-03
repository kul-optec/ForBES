% Broyden method

function [dir, tau0, cache] = Direction_broyden(prob, opt, it, restart, sk, yk, v, cache)

sk = sk(:);
yk = yk(:);

[m, n] = size(v);
v = v(:);

if it == 1 || restart
    dir = -v;
    cache.H = eye(prod(prob.n));
else
    H = cache.H;
    sts = sk'*sk;
    switch opt.modBroyden
    case 3 % absolute value of determinant (guarantees nonsingularity)
        prev_v = cache.prev_v;
        prev_tau = norm(sk)/norm(cache.prev_dir);
        lam = sk'*(H*yk)/sts;
        if abs(lam) < opt.thetaBar
            theta = (1-sign0(lam)*opt.thetaBar)/(1-lam);
            yk = theta*yk - (1-theta)*prev_tau*prev_v;
        end
    otherwise
        error('not implemented');
    end
    Hy = H*yk;
    H = H + (sk-Hy)*(sk'*H)/(sk'*Hy);
    dir = -H*v;
    cache.H = H;
end

cache.prev_v = v;
cache.prev_dir = dir;
tau0 = 1.0;
dir = reshape(dir, m, n);

end
