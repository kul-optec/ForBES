function [cachet, ops] = Cache_SegmentFBE(cache, tau)

if cache.flagLineSearch1 == 0 || cache.flagLineSearch2 == 0
    error('line search data was not precomputed');
end

ops = Ops_Init();

prob = cache.prob;
gam = cache.gam;

cachet = Cache_Init(prob, cache.x + tau*cache.dir1 + (1-tau)*cache.dir2, gam);

fxt = 0;
gradfxt = 0;
if prob.istheref1
    cachet.res1x = cache.res1x + tau*cache.C1dir1 + (1-tau)*cache.C1dir2;
    cachet.gradf1res1x = cache.gradf1res1x + tau*cache.QC1dir1 + (1-tau)*cache.QC1dir2;
    cachet.gradf1x = cache.gradf1x + tau*cache.C1tQC1dir1 + (1-tau)*cache.C1tQC1dir2;
    cachet.f1x = cache.f1x + tau*cache.f1linear + (0.5*tau^2)*cache.f1quad;
    fxt = fxt + cachet.f1x;
    gradfxt = gradfxt + cachet.gradf1x;
end
if prob.istheref2
    cachet.res2x = cache.res2x + tau*cache.C2dir1 + (1-tau)*cache.C2dir2;
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
    cachet.flinx = cache.flinx + tau*cache.lindir1 + (1-tau)*cache.lindir2;
    fxt = fxt + cachet.flinx;
    gradfxt = gradfxt + prob.l;
end
% compute proximal gradient step
cachet.fx = fxt;
cachet.gradfx = gradfxt;
cachet.y = cachet.x - gam*gradfxt;

cachet.flagGradStep = 1;

[cachet.z, cachet.gz] = prob.callg(cachet.y, gam);
ops.proxg = ops.proxg + 1;
ops.g = ops.g + 1;
cachet.FPR = cachet.x-cachet.z;

cachet.flagProxGradStep = 1;

sqnormFPRt = cachet.FPR(:)'*cachet.FPR(:);
cachet.normFPR = sqrt(sqnormFPRt);
cachet.FBE = cachet.fx + cachet.gz - cachet.gradfx(:)'*cachet.FPR(:) + (0.5/gam)*sqnormFPRt;

cachet.flagFBE = 1;

end
