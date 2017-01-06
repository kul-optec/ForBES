% CG-descent

function [dir, cache] = Direction_cgdesc(prob, opt, it, restart, sk, yk, v, cache)

sk = sk(:);
yk = yk(:);

[m, n] = size(v);
v = v(:);

if it == 1 || restart
    dir = -cache_current.gradFBE; % Initially use steepest descent direction
else
    yy = cache_current.gradFBE-cache_previous.gradFBE;
    dy = dir'*yy;
    lambda = 1; %Hager-Zhang proposed lambda = 2 but Kou, Dai found that lambda = 1 is more efficient
    %                 lambda = 2-(dir'*yy)/((dir'*dir)*(yy'*yy));
    beta = ((yy-lambda*dir*(yy'*yy)/dy)'*cache_current.gradFBE)/dy;
    etak = -1/(norm(dir)*min(0.01,norm(cache_current.gradFBE)));
    beta = max(beta,etak);
    dir = -cache_current.gradFBE + beta*dir;
    if dir'*cache_current.gradFBE >= 0 % restart if not descent direction
        dir = -cache_current.gradFBE;
        cntSkip = cntSkip+1;
    end
end

dir = reshape(dir, m, n);

end
