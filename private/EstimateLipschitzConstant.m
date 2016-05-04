% Copyright (C) 2015, Lorenzo Stella and Panagiotis Patrinos
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

function [L, lowerbound] = EstimateLipschitzConstant(prob, opt)
    L = 0;
    lowerbound = 0;
    if isfield(opt, 'Lf')
        L = opt.Lf;
        return;
    end
    eigsOpt.issym = 1;
    eigsOpt.tol = 1e-3;
    if prob.istheref1
        if prob.isthereC1
            if prob.isC1fun
                if prob.isQfun, funHessian = @(x) prob.C1t(prob.Q(prob.C1(x)));
                else funHessian = @(x) prob.C1t(prob.Q*prob.C1(x)); end
            else
                if prob.isQfun, funHessian = @(x) prob.C1'*(prob.Q(prob.C1*x));
                else funHessian = @(x) prob.C1'*(prob.Q*(prob.C1*x)); end
            end
        else
            if prob.isQfun, funHessian = @(x) prob.Q(x);
            else funHessian = @(x) prob.Q*x; end
        end
    end
    if prob.istheref1
        if opt.adaptive
            % compute lower bound for Lipschitz constant
            delta = max(1e-12, prob.x0*1e-6);
            Hx0 = funHessian(prob.x0);
            Hx1 = funHessian(prob.x0+delta);
            L = L + norm(Hx0-Hx1)/norm(delta);
            lowerbound = 1;
        else
            L = L + eigs(funHessian, prob.n, 1, 'LM', eigsOpt);
        end
    end
    if prob.istheref2
        if opt.adaptive || ~isfield(prob, 'Lf2')
            delta = max(1e-12, prob.x0*1e-6);
            if prob.isthereC2
                if prob.isC2fun, C2x0 = prob.C2(prob.x0); else C2x0 = prob.C2*prob.x0; end
                if prob.isC2fun, C2x1 = prob.C2(prob.x0+delta); else C2x1 = prob.C2*(prob.x0+delta); end
            else
                C2x0 = prob.x0;
                C2x1 = prob.x0+delta;
            end
            [~, gradf0] = prob.callf2(C2x0);
            [~, gradf1] = prob.callf2(C2x1);
            L = L + norm(gradf0-gradf1)/norm(delta);
            lowerbound = 1;
        else
            if prob.isthereC2
                if prob.isC2fun
                    funC2tC2 = @(x) prob.C2t(prob.C2(x));
                else
                    funC2tC2 = @(x) prob.C2'*(prob.C2*x);
                end
                sqnormC2 = eigs(funC2tC2, prob.n, 1, 'LM', eigsOpt);
            else
                sqnormC2 = 1;
            end
            L = L + prob.Lf2*sqnormC2;
        end
    end
end
