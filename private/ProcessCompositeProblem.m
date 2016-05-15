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

function prob = ProcessCompositeProblem(prob, opt)
    if ~isfield(prob, 'x0'), error('the starting point x0 must be specified'); end
    prob.x0 = prob.x0;
    prob.n = length(prob.x0);
    if ~any(isfield(prob, {'f1', 'f2'}))
        error('missing f1 and f2');
    end
    if isfield(prob, 'f1')
        if ~isfield(prob.f1, 'isQuadratic') || ~prob.f1.isQuadratic
            error('function f1 must be quadratic');
        end
        if ~isfield(prob.f1, 'makef'), error('value/gradient of f1 is not defined (there is no makef)'); end
        prob.callf1 = prob.f1.makef();
        prob.istheref1 = true;
        if isfield(prob, 'C1')
            if isa(prob.C1, 'function_handle')
                prob.m1 = length(prob.C1(prob.x0));
                if ~isfield(prob, 'C1t') || ~isa(prob.C1t, 'function_handle')
                    error('must specify both C1 and C1t as function handles');
                end
                prob.isC1fun = true;
            else
                prob.m1 = size(prob.C1*prob.x0, 1);
                prob.isC1fun = false;
            end
            prob.isthereC1 = true;
        else
            prob.m1 = prob.n;
            prob.isthereC1 = false;
        end
        [~, ~, prob.Q] = prob.callf1(zeros(prob.m1, 1));
        if isa(prob.Q, 'function_handle'), prob.isQfun = true;
        else prob.isQfun = false; end
        if ~isfield(prob, 'd1'), prob.d1 = zeros(prob.m1, 1); end
        if ~isfield(prob.f1, 'q'), prob.q = zeros(prob.m1, 1);
        else prob.q = prob.f1.q; end
    else
        prob.istheref1 = false;
        prob.isthereC1 = false;
        prob.isC1fun = false;
        prob.isQfun = false;
    end
    if isfield(prob, 'f2')
        if isfield(prob.f2, 'isQuadratic') && prob.f2.isQuadratic
            error('you should provide f2 as f1, since it is quadratic');
        end
        if ~isfield(prob.f2, 'makef'), error('value/gradient of f2 is not defined (there is no makef)'); end
        prob.callf2 = prob.f2.makef();
        prob.istheref2 = true;
        if isfield(prob, 'C2')
            if isa(prob.C2, 'function_handle')
                prob.m2 = length(prob.C2(prob.x0));
                if ~isfield(prob, 'C2t') || ~isa(prob.C2t, 'function_handle')
                    error('must specify both C2 and C2t as function handles');
                end
                prob.isC2fun = true;
            else
                prob.m2 = size(prob.C2, 1);
                prob.isC2fun = false;
            end
            prob.isthereC2 = true;
        else
            prob.m2 = prob.n;
            prob.isthereC2 = false;
        end
        if ~isfield(prob, 'd2'), prob.d2 = zeros(prob.m2, 1); end
        if isfield(prob.f2, 'L'), prob.Lf2 = prob.f2.L; end
        if isfield(prob.f2, 'hasHessian') && prob.f2.hasHessian && opt.useHessian, prob.useHessian = 1;
        else prob.useHessian = 0; end
    else
        prob.istheref2 = false;
        prob.isthereC2 = false;
        prob.isC2fun = false;
    end
    if isfield(prob, 'l')
        prob.istherelin = true;
    else
        prob.istherelin = false;
    end
    if ~isfield(prob, 'g'), error('missing g'); end
    if ~isfield(prob.g, 'makeprox'), error('the prox for the term g you specified is not available'); end
    prob.callg = prob.g.makeprox();
    [prob.Lf, prob.unknownLf] = EstimateLipschitzConstant(prob, opt);
    prob.muf = 0;
    prob.processed = true;
end
