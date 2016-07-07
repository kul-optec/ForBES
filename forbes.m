%FORBES Solver for nonsmooth convex optimization problems.
%
%   Composite problems
%   ------------------
%
%   (1)    minimize f(Cx + d) + g(x)
%
%   We assume that f has Lipschitz continuous gradient,
%   and that g is closed and proper. C is a linear mapping of the
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
%   We assume that f is strongly convex, and that g is closed and proper. 
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
%   on any of the functions.
%
%   Options
%   -------
%
%   In opt the user can specify the behaviour of the algorithm to be used.
%   The following options can be set:
%
%       opt.tol: Tolerance on the optimality condition.
% 
%       opt.maxit: Maximum number of iterations.
%
%		opt.solver: Internal solver to use. Can select between:
%			* 'minfbe' (only for problems where g is convex)
%			* 'zerofpr' (can handle also nonconvex g)
% 
%       opt.method: Algorithm to use. Can select between:
%			* 'bfgs' (BFGS quasi-Newton method)
%           * 'lbfgs' (limited memory BFGS, default).
% 
%       opt.linesearch: Line search strategy to use. Can select between:
%           * 'backtracking' (simple backtracking),
%           * 'backtracking-armijo' (backtracking satisfying Armijo condition),
%           * 'backtracking-nm' (nonmonotone backtracking),
%           * 'lemarechal' (line search for the Wolfe conditions).
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

    if nargin < 3, error('you must provide at least 3 arguments'); end
    if nargin < 4, aff = []; end
    if nargin < 5, constr = []; end
    if nargin < 6, opt = []; end

    prob = MakeProblem(fs, gs, init, aff, constr);
    opt = ProcessOptions(opt);

    switch prob.id

        case 1
            prob = ProcessCompositeProblem(prob, opt);
            preprocess = toc(t0);
            switch opt.solver
                case 'fbs'
                    out = fbs(prob, opt);
                case 'minfbe'
                    out = minfbe(prob, opt);
                case 'zerofpr'
                    out = zerofpr(prob, opt);
            end

        case 2
            [prob, dualprob] = ProcessSeparableProblem(prob, opt);
            preprocess = toc(t0);
            switch opt.solver
                case 'fbs'
                    dualout = fbs(dualprob, opt);
                case 'minfbe'
                    dualout = minfbe(dualprob, opt);
                case 'zerofpr'
                    dualout = zerofpr(dualprob, opt);
            end
            out = GetPrimalOutput(prob, dualprob, dualout);
    end

    out.preprocess = preprocess;
end