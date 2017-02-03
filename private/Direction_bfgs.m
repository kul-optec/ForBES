% BFGS

function [dir, tau0, cache] = Direction_bfgs(prob, opt, it, restart, sk, yk, v, cache)

sk = sk(:);
yk = yk(:);

[m, n] = size(v);
v = v(:);

if it == 1 || restart
    dir = -v;
    cache.H = eye(prod(prob.n));
else
    H = cache.H;
    YSk = yk'*sk;
    Bs = H'*(H*sk);
    sBs = sk'*Bs;
    if YSk > 0
        H = cholupdate(cholupdate(H,yk/sqrt(YSk)),Bs/sqrt(sBs),'-');
    else
        cache.cntSkip = cache.cntSkip+1;
    end
    dir = -linsolve(H,linsolve(H,v,opt.optsL),opt.optsU);
    cache.H = H;
end

tau0 = 1.0;
dir = reshape(dir, m, n);

end
