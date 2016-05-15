function [cache, ops] = CacheGradFBE(cache, gam)

if nargin < 2
    gam = cache.gam;
end

gam0 = cache.gam;

if cache.flagProxGradStep == 0 || gam0 ~= gam
    [cache, ops] = CacheProxGradStep(cache, gam);
else
    ops = OpsInit();
end

if cache.flagGradFBE == 0 || gam0 ~= gam
    prob = cache.prob;
    HFPR = 0;
    if prob.istheref1
        if prob.isthereC1
            if prob.isC1fun, C1FPR = prob.C1(cache.FPR);
            else C1FPR = prob.C1*cache.FPR; end
            if prob.isQfun, Hessf1C1FPR = prob.Q(C1FPR);
            else Hessf1C1FPR = prob.Q*C1FPR; end
            if prob.isC1fun, C1tHessf1C1diff = prob.C1t(Hessf1C1FPR);
            else C1tHessf1C1diff = prob.C1'*Hessf1C1FPR; end
            ops.C1 = ops.C1 + 2;
        else
            if prob.isQfun, C1tHessf1C1diff = prob.Q(cache.FPR);
            else C1tHessf1C1diff = prob.Q*cache.FPR; end
        end
        ops.gradf1 = ops.gradf1 + 1;
        HFPR = HFPR + C1tHessf1C1diff;
    end
    if prob.istheref2
        if prob.isthereC2
            if prob.isC2fun, C2FPR = prob.C2(cache.FPR);
            else C2FPR = prob.C2*cache.FPR; end
            ops.C2 = ops.C2 + 1;
        else
            C2FPR = cache.FPR;
        end
        if prob.useHessian
            HC2FPR = cache.Hessf2res2x*C2FPR;
        else
            % imaginary trick
            res2xepsFPR = cache.res2x + 1e-100i*C2FPR;
            [~, gradf2res2xepsd] = prob.callf2(res2xepsFPR);
            ops.gradf2 = ops.gradf2 + 1;
            HC2FPR = imag(gradf2res2xepsd)/1e-100;
            % forward differences
    %             res2xepsdiff = cache.res2x + 1e-8*C2diff;
    %             [~, gradf2res2xepsd] = prob.callf2(res2xepsdiff);
    %             cnt(4) = cnt(4)+1;
    %             HC2diff = (gradf2res2xepsd-cache.gradf2res2x)/1e-8;
        end
        if prob.isthereC2
            if prob.isC2fun, HFPR = HFPR + prob.C2t(HC2FPR);
            else HFPR = HFPR + (prob.C2'*HC2FPR); end
            ops.C2 = ops.C2 + 1;
        else
            HFPR = HFPR + HC2FPR;
        end
    end
    cache.gradFBE = cache.FPR/gam - HFPR;
    cache.gam = gam;
    cache.flagGradFBE = 1;
end
