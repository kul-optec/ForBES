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

function [prob, dualprob] = ProcessSeparableProblem(prob, opt)
    if isfield(prob, 'b')
        if norm(prob.b) > 0,
            dualprob.l = prob.b;
            dualprob.istherelin = true;
        else
            dualprob.istherelin = false;
        end
        n = length(prob.b);
        if isfield(prob, 'y0'), dualprob.x0 = prob.y0;
        else dualprob.x0 = zeros(n, 1); end
        dualprob.n = n;
    else
        error('you must specify the right hand side b of the equality constraint');
    end
    if ~any(isfield(prob, {'f1', 'f2'}))
        error('missing f1 and f2');
    end
    if isfield(prob, 'f1')
        if ~isfield(prob.f1, 'isConjQuadratic') || ~prob.f1.isConjQuadratic
            error('the conjugate function f1 must be quadratic');
        end
        dualprob.istheref1 = true;
        if ~isfield(prob.f1, 'makefconj'), error('conjugate function of f1 is not defined'); end
        dualprob.isQfun = true;
        if isfield(prob, 'A1')
            dualprob.isthereC1 = true;
            if isa(prob.A1, 'function_handle')
                if ~isfield(prob, 'A1t') || ~isa(prob.A1t, 'function_handle')
                    error('you must specify both A1 and A1t as function handles');
                end
                dualprob.isC1fun = true;
                [dualprob.Q, dualprob.q] = make_quad_conj(prob.f1, prob.A1t(zeros(n, 1)));
                dualprob.C1 = @(x) -prob.A1t(x);
                dualprob.C1t = @(y) -prob.A1(y);
                dualprob.m1 = length(dualprob.C1(dualprob.x0));
            else
                dualprob.isC1fun = false;
                [dualprob.Q, dualprob.q] = make_quad_conj(prob.f1, prob.A1'*zeros(n, 1));
                dualprob.C1 = -prob.A1';
                dualprob.m1 = size(dualprob.C1, 1);
            end
            dualprob.d1 = zeros(dualprob.m1, 1);
        else
            error('you must specify matrix A1 in the constraint');
        end
    else
        dualprob.istheref1 = false;
    end
    if isfield(prob, 'f2')
        if isfield(prob.f2, 'isConjQuadratic') && prob.f2.isConjQuadratic
            error('you should provide f2 as f1, since its conjugate is quadratic');
        end
        dualprob.istheref2 = true;
        if ~isfield(prob.f2, 'makefconj'), error('conjugate function of f2 is not defined'); end
        dualprob.callf2 = prob.f2.makefconj();
        if isfield(prob, 'A2')
            dualprob.isthereC2 = true;
            if isa(prob.A2, 'function_handle')
                dualprob.isthereC2 = true;
                if ~isfield(prob, 'A2t') || ~isa(prob.A2T, 'function_handle')
                    error('you must specify both A2 and A2t as function handles');
                end
                dualprob.isC2fun = true;
                dualprob.C2 = @(x) -prob.A2t(x);
                dualprob.C2t = @(y) -prob.A2(y);
                dualprob.m2 = length(dualprob.C2(dualprob.x0));
            else
                dualprob.isC2fun = false;
                dualprob.C2 = -prob.A2';
                dualprob.m2 = size(dualprob.C2, 1);
            end
            dualprob.d2 = zeros(dualprob.m2, 1);
        else
            error('yout must specify matrix A2 in the constraint');
        end
        if isfield(prob.f2, 'mu'), dualprob.Lf2 = 1/prob.f2.mu; end
        dualprob.useHessian = 0;
    else
        dualprob.istheref2 = false;
    end
    if isfield(prob, 'B')
        if isa(prob.B, 'function_handle'), error('you must specify matrix B as a matrix in the constraint'); end
    else
        error('you must specify matrix B in the constraint');
    end
    mus = sum((prob.B).*(prob.B), 1);
    if (max(mus)-min(mus))/max(mus) > 10*eps, error('B''*B must be a multiple of the identity'); end
    prob.muB = mus(1);
    if ~isfield(prob, 'g'), error('you must specify term g'); end
    if ~isfield(prob.g, 'makeprox'), error('the prox for the term g you specified is not available'); end
    dualprob.callg = make_prox_conj(prob.g, prob.B, prob.muB);
    [dualprob.Lf, dualprob.unknownLf] = EstimateLipschitzConstant(dualprob, opt);
    dualprob.muf = 0;
    dualprob.processed = true;
end

function op = make_prox_conj(g, B, mu)
    prox = g.makeprox();
    op = @(y, gam) call_prox_conj(y, gam, prox, B, mu);
end

function [proxpoint, proxval] = call_prox_conj(y, gam, prox, B, mu)
    mugam = mu*gam;
    [z, v] = prox(-(B'*y)/mugam, 1/mugam);
    Bz = B*z;
    proxpoint = y+gam*Bz;
    proxval = -proxpoint'*Bz - v;
end

function [Q, q] = make_quad_conj(obj, zero)
    fconj = obj.makefconj();
    [~, q] = fconj(zero);
    Q = @(x) hessvec_quad_conj(x, fconj, q);
end
 
function y = hessvec_quad_conj(x, fconj, gradfconj0)
    [~, gradfconjx] = fconj(x);
    y = gradfconjx-gradfconj0;
end
