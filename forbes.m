% FORBES Solvers for nonsmooth, nonconvex optimization problems.
%
%   Composite problems
%   ------------------
%
%   (1)    minimize f(Cx + d) + g(x)
%
%   We assume that f is continuously differentiable, and that g is closed
%   and proper. C is a linear mapping and can be a MATLAB matrix, or any
%   other matrix-like object, which essentially supports matrix-vector
%   products, transposition and 'size'. For example, the SPOT toolbox
%   (http://www.cs.ubc.ca/labs/scl/spot/) can be used to define C.
%
%   out = FORBES(f, g, init, aff, [], sol) solves the problem with the
%   specified f and g. init is the initial value for x, aff is a cell array
%   containing {C, d} (in this order). Last argument sol (optional) is a solver
%   object (more on this later).
%
%   Separable problems
%   ------------------
%
%   (2)    minimize    f(x) + g(z)
%          subject to  Ax + Bz = b
%
%   We assume that f is strongly convex, and that g is closed and proper. A
%   and B are matrix-like objects (just like C in the composite case) and
%   such that B*B' = a*Id for some a > 0.
%
%   out = FORBES(f, g, init, [], constr, opt) solves the specified problem.
%   init is the initial *dual* variable, constr is a cell array defining
%   the constraint, i.e., constr = {A, B, b}. Last argument sol (optional) is a
%   solver object (more on this later).
%
%   General forms of problems
%   -------------------------
%
%   More general forms of problems are formulated as follows, when multiple
%   variables, terms, constraints are involved:
%
%   (1b)     minimize  f_1(C_11 x_1 + ... + C_1M x_M + d_1)
%                    + ...
%                    + f_N(C_N1 x_1 + ... + C_NM x_M + d_N)
%                    + g_M(x_1) + ... + g_M(x_M)
%
%   solved with FORBES(fs, gs, aff) where
%
%         fs = {f_1, ..., f_N}
%         gs = {g_1, ..., g_N}
%         aff = {C_11, ..., C_1M, d_1; ...; C_N1, ..., C_NM, d_N}
%
%   (2b)     minimize   f_1(x_1) + ... + f_N(x_N)
%                     + g_1(z_1) + ... + g_M(z_M)
%          subject to   A_11 x_1 + ... + A_1N x_N + B_1 z_1 = b_1
%                                 [...]
%                       A_M1 x_1 + ... + A_MN x_N + B_M z_M = b_M
%
%   solved with FORBES(fs, gs, [], coeff) where
%
%         fs = {f_1, ..., f_N}
%         gs = {g_1, ..., g_N}
%         coeff = {A_11, ..., A_1N, B_1, b_1; ...; A_M1, ..., A_MN, B_M, b_M}
%
%   Functions and linear mappings
%   -----------------------------
%
%   TODO
%
%   Solvers
%   -------
%
%   TODO
%
%   References
%   ----------
%
%   TODO
%
% Authors: Lorenzo Stella, Panagiotis Patrinos

function out = forbes(fs, gs, init, aff, constr, solver)

    % Fill-in defaults

    if nargin < 4, aff = {}; end
    if nargin < 5, constr = {}; end
    if nargin < 6, solver = forbes.solvers.NAMA(); end

    % Convert single functions to cell arrays

    if ~iscell(fs), fs = {fs}; end
    if ~iscell(gs), gs = {gs}; end

    % Detect problem type (primal or dual)

    ptype = 0;
    if ~isempty(aff)
        ptype = ptype + 1;
    end
    if ~isempty(constr)
        ptype = ptype + 2;
    end

    % Prepare problem

    switch ptype

        case 0

            error('must have either aff or constr');

        case 1 % Solve primal

            n_vars = size(aff, 2)-1;
            [idx_q, idx_n] = split_quadratic(fs);

            [f_q, C_q] = aggregate_terms(fs, aff, idx_q);
            [f_n, C_n] = aggregate_terms(fs, aff, idx_n);

            dims_g = {};
            for i = 1:n_vars
                dims_g{end+1} = size(aff{1, i}, 2);
            end
            if n_vars > 1
                g = forbes.functions.SeparableSum(gs, dims_g);
            elseif n_vars == 1
                g = gs{1};
            else
                error('must have at least one nonsmooth term g');
            end

        case 2 % Solve dual

            error('dual solvers not implemented');

        case 3

            error('cannot have both aff and constr');

        otherwise

            error('unknown problem type');

    end

    % Call solver

    solver.run(f_q, C_q, f_n, C_n, g, init);

    % TODO: should we produce the output differently? probably yes
    % For example, in case we solve the dual, we want the primal solution

    out = solver;

end

function [idx_q, idx_nq] = split_quadratic(fs)
    idx_q = [];
    idx_nq = [];
    for i = 1:length(fs)
        if fs{i}.is_quadratic()
            idx_q = [idx_q, i];
        else
            idx_nq = [idx_nq, i];
        end
    end
end

function [f, C] = aggregate_terms(fs, aff, idx)
    dims = {};
    C = [];
    d = [];
    for i = idx
        dims{end+1} = size(aff{i, 1}, 1);
        C = [C; [aff{i, 1:end-1}]];
        d = [d; aff{i, end}];
    end
    if length(idx) > 1
        s = forbes.functions.SeparableSum(fs(idx), dims);
    elseif length(idx) == 1
        s = fs{idx};
    else
        C = 1.0;
        s = forbes.functions.Zero();
    end
    if norm(d, 'fro') > 0
        f = forbes.functions.ScaleTranslate(s, 1.0, d);
    else
        f = s;
    end
end
