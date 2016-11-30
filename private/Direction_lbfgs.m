% L-BFGS

function [dir, tau0, cache] = Direction_lbfgs(prob, opt, it, restart, sk, yk, v, cache)

sk = full(sk(:));
yk = full(yk(:));
v = full(v(:));

if it == 1 || restart
    dir = -v; % use steepest descent direction initially
    cache.LBFGS_col = 0; % last column of Sk, Yk that was filled in
    cache.LBFGS_mem = 0; % current memory of the method
else
    YSk = yk'*sk;
    if YSk > 0.0
        cache.LBFGS_col = 1+mod(cache.LBFGS_col, opt.memory);
        cache.LBFGS_mem = min(cache.LBFGS_mem+1, opt.memory);
        cache.S(:,cache.LBFGS_col) = sk;
        cache.Y(:,cache.LBFGS_col) = yk;
        cache.YS(cache.LBFGS_col) = YSk;
    else
        cache.cntSkip = cache.cntSkip+1;
    end
    if cache.LBFGS_mem > 0
        H = cache.YS(cache.LBFGS_col)/...
            (cache.Y(:,cache.LBFGS_col)'*cache.Y(:,cache.LBFGS_col));
        dir = lbfgs(cache.S, cache.Y, cache.YS, H, ...
            -v, int32(cache.LBFGS_col), int32(cache.LBFGS_mem));
    else
        dir = -v;
    end
end

tau0 = 1.0;

end
