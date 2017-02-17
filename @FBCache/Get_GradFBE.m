% GET_GRADFBE(cache)
%
%   Returns the gradient of the FBE, the Forward-Backward Envelope, at the
%   point x to which cache refers.

function gradFBE = Get_GradFBE(cache)

if cache.flagGradFBE == true
    gradFBE = cache.gradFBE;
    return;
end

if cache.flagProxGradStep == false
    cache.Get_ProxGradStep();
end

prob = cache.prob;
gam = cache.gam;

HFPR = 0;
if prob.istheref1
    if prob.isthereC1
        C1FPR = prob.C1*cache.FPR;
        QC1FPR = prob.Q(C1FPR);
        C1tQC1diff = prob.C1'*QC1FPR;
        if cache.flagOps
            cache.ops.addC1();
            cache.ops.addC1(); % yes, twice
        end
    else
        C1tQC1diff = prob.Q(cache.FPR);
    end
    if cache.flagOps, cache.ops.addgradf1(); end
    HFPR = HFPR + C1tQC1diff;
end
if prob.istheref2
    if prob.isthereC2
        C2FPR = prob.C2*cache.FPR;
        if cache.flagOps, cache.ops.addC2(); end
    else
        C2FPR = cache.FPR;
    end
    if prob.useHessian
        HC2FPR = cache.Hessf2res2x*C2FPR;
    else
        % imaginary trick
        res2xepsFPR = cache.res2x + 1e-100i*C2FPR;
        [~, gradf2res2xepsd] = prob.callf2(res2xepsFPR);
        if cache.flagOps, cache.ops.addgradf2(); end
        HC2FPR = imag(gradf2res2xepsd)/1e-100;
        % forward differences
%             res2xepsdiff = cache.res2x + 1e-8*C2diff;
%             [~, gradf2res2xepsd] = prob.callf2(res2xepsdiff);
%             cnt(4) = cnt(4)+1;
%             HC2diff = (gradf2res2xepsd-cache.gradf2res2x)/1e-8;
    end
    if prob.isthereC2
        HFPR = HFPR + (prob.C2'*HC2FPR);
        if cache.flagOps, cache.ops.addC2(); end
    else
        HFPR = HFPR + HC2FPR;
    end
end
cache.gradFBE = cache.FPR/gam - HFPR;

cache.flagGradFBE = true;
gradFBE = cache.gradFBE;
