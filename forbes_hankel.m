% FORBES_HANKEL
%
%   FORBES_HANKEL(U, bx, r, mu, opt) solves
%
%       minimize (1/2)||x-b||^2 + mu*sum(svd(H(x)U)
%
%   or equivalently
%
%       minimize (1/2)||x-b||^2 + mu*sum(svd(z))
%       subject to H(x)U - z = 0
%
%   Here, H(x) is a block Hankel matrix formed
%   from x, and b and U come from measurements:
%
%       U is a matrix such that U'U=I
%       x is p by N+1 matrix
%       U is N-r+1 by n matrix
%       H(x) is p(r+1) by N-r+1 Hankel matrix
%

function out = forbes_hankel(U, b, r, mu, tol, freq, opt, out1)

    p = size(b,1);
    N = size(b,2)-1;
    m = (r+1)*p;
    n = size(U,2);

    f = quadLoss(1, b);
    g = nuclearNorm(m, n, mu);
    A = opHankel(N, r, p, U);
    
    if nargin < 8 || isempty(out1)
        y0 = zeros(m,n);
    else
        y0 = out1.y;
    end

    opt.Lf = min(r+1,N-r+1); % squared norm of A
    opt.term = @(prob, it, gam, cache_0, cache_x, ops) terminate_hankel(it, freq, cache_x, b, mu, tol);

    out = forbes(f, g, y0, [], {A, -1, zeros(m, n)}, opt);
    
    G = -A*out.x1;
    nmEY = norm(out.x1-b,'fro')^2/2;
    out.pobj = nmEY + mu*sum(svd(G));
    
end

function flag = terminate_hankel(it, freq, cache, b, mu, tol)
    persistent pobj_best dobj_best;

    if it == 1
        pobj_best = +inf;
        dobj_best = -inf;
    end

    flag = false;

    if it <= 1 || mod(it, freq) ~= 0
        return;
    end

    x = cache.gradf1res1x;
    G = cache.gradf1x;
    nmEY = norm(x-b,'fro')^2/2;
    pobj = nmEY + mu*sum(svd(G));
    if pobj < pobj_best
        pobj_best = pobj;
    end
    dobj = -cache.FBE;
    if dobj > dobj_best
        dobj_best = dobj;
    end
    gap_rel = abs(pobj_best-dobj_best)/max(abs(dobj_best),1);

    if gap_rel < tol
        flag = true;
    end

end
