function [cachet, ops] = Cache_LineFBE(cache, tau, mode, cachet)
% Computes the 'directional' FBE, i.e., FBE(x+tau*d) and its derivative with
% respect to tau, if requested. Here x = cache.x, d = cache.dir.
% If cachet (4th argument) is provided, then skips precomputing the data
% that has already been stored in cachet.
%
% If mode == 1, then compute only FBE(x+tau*d) and put it into cachet.
% If mode == 2, compute only dFBE(x+tau*d), the directional derivative.
% If mode == 3, compute both FBE and dFBE at x+tau*d.

if cache.flagLineSearch1 ~= 1
    error('line search data was not precomputed');
end

ops = Ops_Init();

prob = cache.prob;
gam = cache.gam;

if nargin < 4
    cachet = Cache_Init(prob, cache.x + tau*cache.dir1, gam);
end

if nargin < 4 || cachet.flagGradStep == 0
    fxt = 0;
    gradfxt = 0;
    if prob.istheref1
        cachet.res1x = cache.res1x + tau*cache.C1dir1;
        cachet.gradf1res1x = cache.gradf1res1x + tau*cache.QC1dir1;
        cachet.gradf1x = cache.gradf1x + tau*cache.C1tQC1dir1;
        cachet.f1x = cache.f1x + tau*cache.f1linear1 + (0.5*tau^2)*cache.f1quad1;
        fxt = fxt + cachet.f1x;
        gradfxt = gradfxt + cachet.gradf1x;
    end
    if prob.istheref2
        cachet.res2x = cache.res2x + tau*cache.C2dir1;
        if prob.useHessian
            [f2xt, gradf2res2xt, cachet.Hessf2res2x] = prob.callf2(cachet.res2x);
        else
            [f2xt, gradf2res2xt] = prob.callf2(cachet.res2x);
            cachet.gradf2res2x = gradf2res2xt;
        end
        ops.f2 = ops.f2 + 1;
        ops.gradf2 = ops.gradf2 + 1;
        if prob.isthereC2
            gradf2xt = prob.C2'*gradf2res2xt;
            ops.C2 = ops.C2 + 1;
        else
            gradf2xt = gradf2res2xt;
        end
        fxt = fxt + f2xt;
        gradfxt = gradfxt + gradf2xt;
    end
    if prob.istherelin
        cachet.flinx = cache.flinx + tau*cache.lindir1;
        fxt = fxt + cachet.flinx;
        gradfxt = gradfxt + prob.l;
    end
    % compute proximal gradient step
    cachet.fx = fxt;
    cachet.gradfx = gradfxt;
    cachet.y = cachet.x - gam*gradfxt;

    cachet.flagGradStep = 1;
end

if nargin < 4 || cachet.flagProxGradStep == 0
    [cachet.z, cachet.gz] = prob.callg(cachet.y, gam);
    ops.proxg = ops.proxg + 1;
    ops.g = ops.g + 1;
    cachet.FPR = cachet.x-cachet.z;

    cachet.flagProxGradStep = 1;
end

if mode == 1 || mode == 3
    sqnormFPRt = cachet.FPR(:)'*cachet.FPR(:);
    cachet.normFPR = sqrt(sqnormFPRt);
    cachet.FBE = cachet.fx + cachet.gz - cachet.gradfx(:)'*cachet.FPR(:) + (0.5/gam)*sqnormFPRt;

    cachet.flagFBE = 1;
end

if mode >= 2
    Hdir1 = 0;
    if prob.istheref1
        Hdir1 = Hdir1 + cache.C1tQC1dir1;
    end
    if prob.istheref2
        if prob.useHessian
            HC2dir1 = cachet.Hessf2res2x*cache.C2dir1;
        else
            res2xtepsdir1 = cachet.res2x + 1e-100i*cache.C2dir1;
            [~, gradf2res2xtepsdir1] = prob.callf2(res2xtepsdir1);
            ops.gradf2 = ops.gradf2 + 1;
            HC2dir1 = imag(gradf2res2xtepsdir1)/1e-100;
        end
        if prob.isthereC2
            Hdir1 = Hdir1 + (prob.C2'*HC2dir1);
            ops.C2 = ops.C2 + 1;
        else
            Hdir1 = Hdir1 + HC2dir1;
        end
    end
    cachet.dFBE = (cachet.FPR(:)'*cache.dir1(:))/gam - cachet.FPR(:)'*Hdir1(:);
end

end
