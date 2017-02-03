% restarted Broyden method

function [dir, tau0, cache] = Direction_rbroyden(prob, opt, it, restart, sk, yk, v, cache)

sk = sk(:);
yk = yk(:);

[m, n] = size(v);
v = v(:);

if it == 1 || restart
    dir = -v;
    cache.S = [];
    cache.HY = [];
    cache.shy = [];
else
    sts = sk'*sk;
    yts = yk'*sk;
    yty = yk'*yk;
    if opt.initialScaling
        h = yts/yty;
    else
        h = 1;
    end
    dir = -h*v;
    hy = h*yk;
    for i=1:size(cache.S, 2)
      	hy = hy + ((hy'*cache.S(:,i))/cache.shy(i)) * (cache.S(:,i)-cache.HY(:,i));
      	dir = dir + ((dir'*cache.S(:,i))/cache.shy(i)) * (cache.S(:,i)-cache.HY(:,i));
    end
    shy = hy'*sk;
    % damping
    switch opt.modBroyden
    case 1 % enforces positive curvature along sk
        % compute theta_{k-1}; if not 1 then update HYk = H_{k-1}\tilde y_{k-1}
        prev_v = cache.prev_v;
        prev_tau = norm(sk)/norm(cache.prev_dir);
        delta = -prev_tau*(sk'*prev_v); % delta = \delta_{k-1} = <B_{k-1}s_{k-1}, s_{k-1}> = -tau*<s_{k-1},Rw^{k-1}>  (Rw^k = Rxold)
        if yts < opt.deltaCurvature*abs(delta)
            theta = (1-sign0(delta)*opt.deltaCurvature)*delta/(delta-yts);
            hy = (1-theta)*sk + theta*hy;
            shy = (1-theta)*sts + theta*shy;
        end
    case 3 % nonsingularity
    		% compute theta_{k-1}; if not 1 then update HYk = H_{k-1}\tilde y_{k-1}
    		gam  = shy/sts;
    		if abs(gam) < opt.thetaBar
    			theta = (1-sign0(gam)*opt.thetaBar)/(1-gam);
    			hy = (1-theta)*sk + theta*hy;
    			shy = (1-theta)*sts + theta*shy;
    		end
    otherwise
        error('not implemented');
    end
    dir = dir + ((dir'*sk)/shy) * (sk-hy);
    % update buffer
    if size(cache.S,2) < opt.memory
      	cache.S = [cache.S, sk];
      	cache.HY = [cache.HY, hy];
      	cache.shy = [cache.shy, shy];
    else
      	cache.S  = [];
      	cache.HY = [];
      	cache.shy = [];
    end
end

cache.prev_v = v;
cache.prev_dir = dir;
tau0 = 1.0;
dir = reshape(dir, m, n);

end
