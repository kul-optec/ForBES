% Copyright (C) 2015-2016, Lorenzo Stella and Panagiotis Patrinos
%
% This file is part of ForBES.
%
% ForBES is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ForBES is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with ForBES. If not, see <http://www.gnu.org/licenses/>.

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
            C1FPR = prob.C1*cache.FPR;
            QC1FPR = prob.Q(C1FPR);
            C1tQC1diff = prob.C1'*QC1FPR;
            ops.C1 = ops.C1 + 2;
        else
            C1tQC1diff = prob.Q(cache.FPR);
        end
        ops.gradf1 = ops.gradf1 + 1;
        HFPR = HFPR + C1tQC1diff;
    end
    if prob.istheref2
        if prob.isthereC2
            C2FPR = prob.C2*cache.FPR;
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
            HFPR = HFPR + (prob.C2'*HC2FPR);
            ops.C2 = ops.C2 + 1;
        else
            HFPR = HFPR + HC2FPR;
        end
    end
    cache.gradFBE = cache.FPR/gam - HFPR;
    cache.gam = gam;
    cache.flagGradFBE = 1;
end
