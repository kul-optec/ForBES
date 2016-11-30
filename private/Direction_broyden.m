% Broyden method

function [dir, tau0, cache] = Direction_broyden(prob, opt, it, restart, sk, yk, v, cache)

sk = sk(:);
yk = yk(:);
v = v(:);

if it == 1 || restart % || mod(it, 5) == 1
    dir = -v;
    cache.R = eye(prod(prob.n));
else
    R = cache.R;
    sts = sk'*sk;
    if opt.bopt == 1 % enforces nonzero determinant
        sig = 0.1;
        prev_v = cache.prev_v;
        prev_tau = cache.prev_tau;
        gammak = ((R*yk)'*sk)/sts;
        if abs(gammak) < sig,
            theta = (1-sgn(gammak)*sig)/(1-gammak);
            yk = theta*yk - (1-theta)*prev_tau*prev_v;
        end
    elseif opt.bopt == 2 % enforces positive curvature along sk
        sig = 0.5;
        prev_v = cache.prev_v;
        prev_tau = cache.prev_tau;
        sty = sk'*yk;
        stv = sk'*prev_v;
        if sty < sig*prev_tau*abs(stv)
            theta = (1+sgn(stv)*sig)*prev_tau*stv/(prev_tau*stv + sty);
            yk = theta*yk - (1-theta)*prev_tau*prev_v;
        end
    end
    Ry = R*yk;
    R = R + (sk-Ry)*(sk'*R)/(sk'*Ry);
    dir = -R*v;
    cache.R = R;
end

tau0 = 1.0;

end
