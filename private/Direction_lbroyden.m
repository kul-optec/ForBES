% limited memory Broyden method

function [dir, tau0, cache] = Direction_lbroyden(prob, opt, it, restart, sk, yk, v, cache)

sk = sk(:);
yk = yk(:);

[m, n] = size(v);
v = v(:);

if it == 1 || restart
    dir = -v; % use steepest descent direction initially
    cache.S = []; % stores vectors sk
    cache.Y = []; % stores vectors yk
    cache.W = []; % stores columns of H0*Y-S
    cache.StY = []; % stores inner products <sk,yk>
    cache.M = []; % these two are m-by-m where m is the memory of the method
    cache.LBroyden_mem = 0;
else
    % damping
    if opt.modBroyden == 2 % enforces positive curvature along sk
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
    wk = delta*yk - sk;
    if cache.LBroyden_mem == opt.memory, idx0 = 2;
    else idx0 = 1; cache.LBroyden_mem = cache.LBroyden_mem+1; end
    S0 = cache.S(:,idx0:end);
    Y0 = cache.Y(:,idx0:end);
    W0 = cache.W(:,idx0:end);
    StY0 = cache.StY(idx0:end,idx0:end);
    M0 = cache.M(idx0:end,idx0:end);
    % update matrices S, Y, W, StY, M
    cache.S = [S0, sk];
    cache.Y = [Y0, yk];
    cache.W = [W0, wk];
    if isempty(Y0), cache.StY = sk'*yk;
    else cache.StY = [[StY0; sk'*Y0], cache.S'*yk]; end
    if isempty(S0), cache.M = 0;
    else cache.M = [[M0; S0(:,end)'*S0], zeros(cache.LBroyden_mem, 1)]; end
    K = delta * cache.StY - cache.M;
    % compute direction
    dir = delta*(cache.W * (K\(cache.S'*v)) - v);
end

tau0 = 1.0;
dir = reshape(dir, m, n);

end
