% BFGS

function [dir, tau0, cache] = Direction_bfgs(prob, opt, it, restart, sk, yk, v, cache)

sk = sk(:);
yk = yk(:);
v = v(:);

if it == 1 || restart
    dir = -v;
    cache.R = eye(prod(prob.n));
else
    R = cache.R;
    YSk = yk'*sk;
    Bs = R'*(R*sk);
    sBs = sk'*Bs;
    if YSk > 0
        R = cholupdate(cholupdate(R,yk/sqrt(YSk)),Bs/sqrt(sBs),'-');
    else
        cache.cntSkip = cache.cntSkip+1;
    end
    dir = -linsolve(R,linsolve(R,v,opt.optsL),opt.optsU);
    cache.R = R;
end

tau0 = 1.0;

end
