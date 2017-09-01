function [cachet, ops] = Get_CacheSegment(cache, tau)

if ~cache.flagLineSearch1 || ~cache.flagLineSearch2
    error('line search data was not precomputed');
end

prob = cache.prob;
gam = cache.gam;

cachet = forbes.fbe.FBCache(cache.prob, cache.x + tau*cache.dir1 + (1-tau)*cache.dir2, cache.gam, cache.ops);

fxt = 0;
gradfxt = 0;
if prob.istheref1
    cachet.res1x = cache.res1x + tau*cache.C1dir1 + (1-tau)*cache.C1dir2;
    cachet.gradf1res1x = cache.gradf1res1x + tau*cache.QC1dir1 + (1-tau)*cache.QC1dir2;
    cachet.gradf1x = cache.gradf1x + tau*cache.C1tQC1dir1 + (1-tau)*cache.C1tQC1dir2;
    cachet.f1x = cache.f1x + tau*cache.f1linear1 + (1-tau)*cache.f1linear2 + (0.5*tau^2)*cache.f1quad1 + (0.5*(1-tau)^2)*cache.f1quad2 + tau*(1-tau)*cache.f1cross;
    fxt = fxt + cachet.f1x;
    gradfxt = gradfxt + cachet.gradf1x;
end
if prob.istheref2
    cachet.res2x = cache.res2x + tau*cache.C2dir1 + (1-tau)*cache.C2dir2;
    if prob.useHessian
        [gradf2res2xt, f2xt, cachet.Hessf2res2x] = prob.f2.gradient(cachet.res2x);
    else
        [gradf2res2xt, f2xt] = prob.f2.gradient(cachet.res2x);
        cachet.gradf2res2x = gradf2res2xt;
    end
    if cache.flagOps
        cache.ops.addf2();
        cache.ops.addgradf2();
    end
    if prob.isthereC2
        gradf2xt = prob.C2'*gradf2res2xt;
        if cache.flagOps, cache.ops.addC2(); end
    else
        gradf2xt = gradf2res2xt;
    end
    fxt = fxt + f2xt;
    gradfxt = gradfxt + gradf2xt;
end
if prob.istherelin
    cachet.flinx = cache.flinx + tau*cache.lindir1 + (1-tau)*cache.lindir2;
    fxt = fxt + cachet.flinx;
    gradfxt = gradfxt + prob.lin;
end
% compute proximal gradient step
cachet.fx = fxt;
cachet.gradfx = gradfxt;
cachet.y = cachet.x - gam*gradfxt;

cachet.flagGradStep = true;

if prob.isthereD
    mugam = prob.mu*gam;
    [z, cachet.gz] = prob.g.prox(prob.D*cachet.y, mugam);
    cachet.z = cachet.y + prob.D'*(z - prob.D*cachet.y)/prob.mu;
else
    [cachet.z, cachet.gz] = prob.g.prox(cachet.y, gam);
end
if cache.flagOps
    cache.ops.addproxg();
    cache.ops.addg();
end
cachet.FPR = cachet.x-cachet.z;

cachet.flagProxGradStep = true;

sqnormFPRt = cachet.FPR(:)'*cachet.FPR(:);
cachet.normFPR = sqrt(sqnormFPRt);
cachet.FBE = cachet.fx + cachet.gz - cachet.gradfx(:)'*cachet.FPR(:) + (0.5/gam)*sqnormFPRt;

cachet.flagFBE = true;
