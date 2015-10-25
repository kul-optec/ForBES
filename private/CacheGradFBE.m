function [cache, ops] = CacheGradFBE(cache, gam)

ops = OpsInit();

if cache.flagProxGradStep == 0 || cache.gam ~= gam
    [cache, ops] = CacheProxGradStep(cache, gam);
end

prob = cache.prob;

Hdiff = 0;
if prob.istheref1
    if prob.isthereC1
        if prob.isC1fun, C1diff = prob.C1(cache.diff);
        else C1diff = prob.C1*cache.diff; end
        if prob.isQfun, QC1diff = prob.Q(C1diff);
        else QC1diff = prob.Q*C1diff; end
        if prob.isC1fun, C1tQC1diff = prob.C1t(QC1diff);
        else C1tQC1diff = prob.C1'*QC1diff; end
        ops.C1 = ops.C1 + 2;
    else
        if prob.isQfun, C1tQC1diff = prob.Q(cache.diff);
        else C1tQC1diff = prob.Q*cache.diff; end
    end
    ops.Q = ops.Q + 2;
    Hdiff = Hdiff + C1tQC1diff;
end
if prob.istheref2
    if prob.isthereC2
        if prob.isC2fun, C2diff = prob.C2(cache.diff);
        else C2diff = prob.C2*cache.diff; end
        ops.C2 = ops.C2 + 1;
    else
        C2diff = cache.diff;
    end
    if prob.useHessian
        HC2diff = cache.Hessf2res2x*C2diff;
    else
        % imaginary trick
        res2xepsdiff = cache.res2x + 1e-100i*C2diff;
        [~, gradf2res2xepsd] = prob.callf2(res2xepsdiff);
        ops.gradf2 = ops.gradf2 + 1;
        HC2diff = imag(gradf2res2xepsd)/1e-100;
        % forward differences
%             res2xepsdiff = cache.res2x + 1e-8*C2diff;
%             [~, gradf2res2xepsd] = prob.callf2(res2xepsdiff);
%             cnt(4) = cnt(4)+1;
%             HC2diff = (gradf2res2xepsd-cache.gradf2res2x)/1e-8;
    end
    if prob.isthereC2
        if prob.isC2fun, Hdiff = Hdiff + prob.C2t(HC2diff);
        else Hdiff = Hdiff + (prob.C2'*HC2diff); end
        ops.C2 = ops.C2 + 2;
    else
        Hdiff = Hdiff + HC2diff;
    end
end

cache.gradFBE = (Hdiff - cache.diff/gam);
cache.gam = gam;

cache.flagGradFBE = 1;
