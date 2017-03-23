% BFGS, naive implementation intended for debugging purposes
% (direct update & backslash)

function [dir, tau0, cache] = Direction_bfgs_naive(prob, opt, it, restart, sk, yk, v, cache)

sk = sk(:);
yk = yk(:);

[m, n] = size(v);
v = v(:);

if it == 1 || restart
    dir = -v;
    cache.B = eye(prod(prob.n));
else
    B = cache.B;
    YSk = yk'*sk;
    if YSk > 0
        Bs = B*sk;
        sBs = sk'*Bs;
        B = B + (yk*yk')/YSk - (Bs*Bs')/sBs;
    else
        cache.cntSkip = cache.cntSkip+1;
    end
    dir = -B\v;
    cache.B = B;
end

tau0 = 1.0;
dir = reshape(dir, m, n);

end
