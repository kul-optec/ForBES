%FORBES Solver for nonsmooth convex optimization problems.
%
%   Composite problems
%   ------------------
%
%   (1)    minimize f(Cx + d) + g(x)
%
%   We assume that f is convex smooth with Lipschitz continuous gradient,
%   and that g is closed, proper, convex. C is a linear mapping of the
%   appropriate dimension.
%
%   out = FORBES(f, g, init, aff, [], opt) solves the problem with the
%   specified f and g. init is the initial value for x, aff is a cell array
%   containing {C, d} (in this order). opt is a structure defining the
%   options for the solver (more on this later).
%
%   Separable problems
%   ------------------
%
%   (2)    minimize    f(x) + g(z)
%          subject to  Ax + Bz = b
%
%   We assume that f is strongly convex, and that g is closed, proper,
%   convex. 
%
%   out = FORBES(f, g, init, [], constr, opt) solves the specified problem.
%   init is the initial *dual* variable, constr is a cell array defining
%   the constraint, i.e., constr = {A, B, b}. the options are specified in
%   the opt structure (more on this later).
%
%   Functions and linear mappings
%   -----------------------------
%
%   Functions f and g in the cost can be selected in a library of functions
%   available in the "library" directory inside of FORBES directory. Linear
%   mappings (C in problem (1) and A, B in problem (2) above) can either be
%   MATLAB's matrices or can themselves be picked from a library of
%   standard linear operators.
%
%   For example, to define f and g:
%
%       f = logLoss() % logistic loss function
%       g = l1Norm() % l1 regularization term
%
%   Consider looking into the "library" directory for specific information
%   any of the functions.
%
%   Options
%   -------
%
%   In opt the user can specify the behaviour of the algorithm to be used.
%   The following options can be set:
%
%       opt.tol: Tolerance on the optimality condition. (Default: 1e-8).
% 
%       opt.maxit: Maximum number of iterations. (Default: 10*n, with n
%           being the number of variables).
% 
%       opt.method: Algorithm to use. Can select between:
%           * 'sd' (steepest descent),
%           * 'lbfgs' (limited memory BFGS, default),
%           * 'cg-desc', 'cg-prp', 'cg-dyhs' (various CG algorithms),
%           * 'bb' (Barzilai-Borwein).
% 
%       opt.variant: Variant of the method to use. Select between:
%           * 'basic', the basic algorithm,
%           * 'global', the global variant (default),
%           * 'fast', the fast variant.
% 
%       opt.linesearch: Line search strategy to use. Can select between:
%           * 'armijo' (default for 'sd'),
%           * 'nonmonotone-armijo' (default for 'bb'),
%           * 'hager-zhang' (default for the rest),
%           * 'lemarechal',
%           * 'fletcher'.
%
% Authors: Lorenzo Stella (lorenzo.stella -at- imtlucca.it)
%          Panagiotis Patrinos (panagiotis.patrinos -at- imtlucca.it)
%
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

function out = forbes(fs, gs, init, aff, constr, opt)
    t0 = tic();
    if nargin < 4, aff = []; end
    if nargin < 5, constr = []; end
    if nargin < 6, opt = []; end
    opt = ProcessOptions(opt);
    prob = MakeProb(fs, gs, init, aff, constr, opt);
    switch prob.id
        case 1
            prob = ProcessCompositeProblem(prob, opt);
            preprocess = toc(t0);
            if opt.method == 0, out = fbs(prob, opt);
            else out = minfbe(prob, opt); end
        case 2
            [prob, dualprob] = ProcessSeparableProblem(prob, opt);
            preprocess = toc(t0);
            if opt.method == 0, dualout = fbs(dualprob, opt);
            else dualout = minfbe(dualprob, opt); end
            out = GetPrimalOutput(prob, dualprob, dualout);
    end
    out.preprocess = preprocess;
end

function prob = MakeProb(fs, gs, init, aff, constr, opt)
    M = length(fs);
    N = length(gs);
    if M > 2
        error('only the sum of two functions is currently supported');
    end
    if ~isa(fs, 'cell'), fs = {fs}; end
    if isempty(constr)
        flagconstr = 0;
    else
        if ~isa(constr, 'cell')
            error('the constraint must be a cell array');
        end
        flagconstr = 1;
        if length(constr) ~= N+M+1
            error('must have as many terms in the constraints as f and g functions, plus the rhs');
        end
    end
    if isempty(aff)
        flagaff = 0;
        flagd = 0;
        ns = length(init);
    else
        if isa(aff, 'double') || isa(aff, 'struct')
            aff = {aff};
        end
        if ~isa(aff, 'cell')
            error('the list of affine maps must be a cell array, a matrix or a structure');
        end
        flagaff = 1;
        if length(aff) == N, flagd = 0;
        elseif length(aff) == N+1, flagd = 1;
        else error('must have as many blocks of variables as g functions, if any');
        end
        ns = zeros(1,N);
        for i=1:N
            if ismatrix(aff{1,i}), ns(i) = size(aff{1,i},2);
            else ns(i) = aff{1,i}.n; end
        end
    end
    prob.id = flagconstr+1;
    switch prob.id
        case 1
            prob.x0 = init;
            for i = 1:M
                if isfield(fs{i}, 'isQuadratic') && fs{i}.isQuadratic
                    prob.f1 = fs{i};
                    if flagaff
                        op = separableSumLinear(aff(i,1:N));
                        if isstruct(op)
                            prob.C1 = op.makeop();
                            prob.C1t = op.makeadj();
                        else
                            prob.C1 = op;
                        end
                        if flagd, prob.d1 = -aff{i,N+1}; end
                    end
                end
                if ~isfield(fs{i}, 'isQuadratic') || ~fs{i}.isQuadratic
                    prob.f2 = fs{i};
                    if flagaff
                        op = separableSumLinear(aff(i,1:N));
                        if isstruct(op)
                            prob.C2 = op.makeop();
                            prob.C2t = op.makeadj();
                        else
                            prob.C2 = op;
                        end
                        if flagd, prob.d2 = -aff{i,N+1}; end
                    end
                end
            end
            if N == 1
                if isa(gs, 'cell'), prob.g = gs{1};
                else prob.g = gs; end
            elseif N > 1
                prob.g = separableSum(gs, ns);
            end
        case 2
            prob.y0 = init;
            if flagaff
                error('cannot have both constraints and affine mappings');
            end
            for i = 1:M
                if isfield(fs{i}, 'isConjQuadratic') && fs{i}.isConjQuadratic
                    prob.f1 = fs{1};
                    if isstruct(constr{i})
                        prob.A1 = constr{i}.makeop();
                        prob.A1t = constr{i}.makeadj();
                    else
                        prob.A1 = constr{i};
                    end
                end
                if ~isfield(fs{i}, 'isConjQuadratic') || ~fs{i}.isConjQuadratic
                    prob.f2 = fs{i};
                    if isstruct(constr{i})
                        prob.A2 = constr{i}.makeop();
                        prob.A2t = constr{i}.makeadj();
                    else
                        prob.A2 = constr{i};
                    end
                end
            end
            if N == 1
                if isa(gs, 'cell'), prob.g = gs{1};
                else prob.g = gs; end
            elseif N > 1
                prob.g = separableSum(gs, ns);
            end
            prob.B = horzcat(constr{M+1:M+N});
            prob.b = constr{M+N+1};
    end
end

function op = separableSumLinear(linops)
    N = length(linops);
    % for a single linear operator, return the operator itself
    if N == 1
        if isa(linops, 'cell') op = linops{1};
        else op = linops; end
        return;
    end
    flagstruct = zeros(1,N);
    n = 0;
    for i = 1:N
        if isstruct(linops{i})
            flagstruct(i) = 1;
            m = linops{i}.m;
            n = n + linops{i}.n;
            linops{i}.op = linops{i}.makeop();
            linops{i}.adj = linops{i}.makeadj();
        else
            m = size(linops{i},1);
            n = n + size(linops{i},2);
        end
    end
    % if all operators are matrices, simply stack them horizontally
    if flagstruct == 0, op = horzcat(linops{:}); return; end
    % otherwise make a structure
    op.m = m;
    op.n = n;
    op.makeop = @() @(x) call_separableSumLinear(x, linops, N, m, flagstruct);
    op.makeadj = @() @(y) call_separableSumLinear_adj(y, linops, N, flagstruct);
end

function y = call_separableSumLinear(x, linops, N, m, flags)
    y = zeros(m,1);
    k = 0;
    for i=1:N
        if flags(i)
            y = y + linops{i}.op(x(k+1:k+linops{i}.n));
        else
            y = y + linops{i}*x(k+1:k+linops{i}.n);
        end
    end
end

function x = call_separableSumLinear_adj(y, linops, N, flags)
    x = [];
    for i=1:N
        if flags(i)
            x = [x; linops{i}.adj(y)];
        else
            x = [x; linops{i}*y];
        end
    end
end
