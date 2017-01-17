function gradFBE = Get_GradFBE(cache, gam)

if nargin < 2
    gam = cache.gam;
end

gam0 = cache.gam;

if ~cache.flagProxGradStep || gam0 ~= gam
    cache.ProxGradStep(gam);
end

if ~cache.flagGradFBE || gam0 ~= gam
    prob = cache.prob;
    HFPR = 0;
    if prob.istheref1
        if prob.isthereC1
            C1FPR = prob.C1*cache.FPR;
            QC1FPR = prob.Q(C1FPR);
            C1tQC1diff = prob.C1'*QC1FPR;
            cache.ops.addC1();
            cache.ops.addC1(); % yes, twice
        else
            C1tQC1diff = prob.Q(cache.FPR);
        end
        cache.ops.addgradf1();
        HFPR = HFPR + C1tQC1diff;
    end
    if prob.istheref2
        if prob.isthereC2
            C2FPR = prob.C2*cache.FPR;
            cache.ops.addC2();
        else
            C2FPR = cache.FPR;
        end
        if prob.useHessian
            HC2FPR = cache.Hessf2res2x*C2FPR;
        else
            % imaginary trick
            res2xepsFPR = cache.res2x + 1e-100i*C2FPR;
            [~, gradf2res2xepsd] = prob.callf2(res2xepsFPR);
            cache.ops.addgradf2();
            HC2FPR = imag(gradf2res2xepsd)/1e-100;
            % forward differences
    %             res2xepsdiff = cache.res2x + 1e-8*C2diff;
    %             [~, gradf2res2xepsd] = prob.callf2(res2xepsdiff);
    %             cnt(4) = cnt(4)+1;
    %             HC2diff = (gradf2res2xepsd-cache.gradf2res2x)/1e-8;
        end
        if prob.isthereC2
            HFPR = HFPR + (prob.C2'*HC2FPR);
            cache.ops.addC2();
        else
            HFPR = HFPR + HC2FPR;
        end
    end
    cache.gradFBE = cache.FPR/gam - HFPR;
    cache.gam = gam;
    cache.flagGradFBE = true;
end

gradFBE = cache.gradFBE;
