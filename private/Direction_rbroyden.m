% restarted Broyden method

function [dir, tau0, cache] = Direction_rbroyden(prob, opt, it, restart, sk, yk, v, cache)

sk = sk(:);
yk = yk(:);
v = v(:);

if it == 1 || restart || mod(it, opt.memory) == 1
    dir = -v;
    cache.LBroyden_mem = 0;
    cache.S = [];
    cache.Y = [];
    cache.W = [];
else
    % damping
    if opt.bopt == 2 % enforces positive curvature along sk
        sig = 0.1;
        prev_v = cache.prev_v;
        prev_tau = cache.prev_tau;
        sty = sk'*yk;
        stv = sk'*prev_v;
        if sty < sig*prev_tau*abs(stv)
            theta = (1+sign(stv)*sig)*prev_tau*stv/(prev_tau*stv + sty);
            yk = theta*yk - (1-theta)*prev_tau*prev_v;
        end
    end
    delta = 1; % diagonal of H0
    if cache.LBroyden_mem == opt.memory, idx0 = 2;
    else idx0 = 1; cache.LBroyden_mem = cache.LBroyden_mem+1; end
    S0 = cache.S(:,idx0:end);
    Y0 = cache.Y(:,idx0:end);
    W0 = cache.W(:,idx0:end);
    w = delta*yk;
    dir = -delta*v;
    for j = 1:size(W0, 2)
        w = w + (S0(:,j)'*w)*W0(:,j);
        dir = dir + (S0(:,j)'*dir)*W0(:,j);
    end
    wk = (sk-w)/(sk'*w);
    dir = dir + (sk'*dir)*wk;
    % update matrices S, Y, W
    cache.S = [S0, sk];
    cache.Y = [Y0, yk];
    cache.W = [W0, wk];
end

tau0 = 1.0;

end
