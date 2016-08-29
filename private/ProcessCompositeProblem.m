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

function prob = ProcessCompositeProblem(prob, opt)
    if ~isfield(prob, 'x0'), error('the starting point x0 must be specified'); end
    prob.n = size(prob.x0);
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
            prob.m1 = [size(prob.C1, 1), 1];
            prob.isthereC1 = true;
        else
            prob.m1 = prob.n;
            prob.isthereC1 = false;
        end
        if ~isfield(prob, 'd1'), prob.d1 = zeros(prob.m1); end
        [~, q] = prob.callf1(zeros(prob.m1));
        prob.Q = @(x) HessianQuadratic(x, prob.callf1, q);
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
            prob.m2 = [size(prob.C2, 1), 1];
            prob.isthereC2 = true;
        else
            prob.m2 = prob.n;
            prob.isthereC2 = false;
        end
        if ~isfield(prob, 'd2'), prob.d2 = zeros(prob.m2); end
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
    if isfield(opt, 'Lf')
        prob.Lf = opt.Lf;
        prob.adaptive = 0;
    else
        prob.Lf = EstimateLipschitzConstant(prob);
        prob.adaptive = 1;
    end
    prob.muf = 0;
    prob.processed = true;
end
