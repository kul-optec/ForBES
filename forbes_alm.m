%FORBES_ALM Solver for nonsmooth optimization problems.
%
%   FORBES_ALM(f, g, h, F, init, opt, inn_init, inn_opt) solves problems of the form
%
%       minimize f(x) + g(x) + h(F(x))
%
%   We assume that f is smooth with Lipschitz continuous gradient, that
%   g is a closed, proper function and that h is a closed, proper, convex
%   function. Both g and h are assumed to have an easily computable proximal
%   mapping. F is a linear mapping.
%
%   Parameter init is the starting dual point, opt is a structure
%   containing options for the augmented Lagrangian method, inn_opt is a
%   structure containing options for the inner solver.
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

function out = forbes_alm(f, g, h, F, init, opt, inn_init, inn_opt)

if nargin < 5 || isempty(init), init = zeros(F.m, 1); end
if nargin < 6, opt = []; end
if nargin < 7 || isempty(inn_init), inn_init = zeros(F.n,1); end
if nargin < 8, inn_opt = []; end

% fill-in missing options with defaults
opt = default_opt(opt);
inn_opt = default_inner_opt(inn_opt);

if opt.display >= 2
    fprintf('%6s%11s%11s%11s%11s\n', 'iter', 'res', 'penalty', 'inner tol', 'inner it');
end

y = init;
x = inn_init;
res = zeros(1,opt.maxit);
callF = F.makeop();

tot_inn_it = 0;
tot_ops = OpsInit();

% Algorithm 17.4 Nocedal
r = 10;
inn_opt.tol = 1/r;
eta = 0.1258925;

inn_linop = stackOp({F, diagOp(F.n)});

if ~isfield(opt, 'sqOpNorm')
    linop_op = inn_linop.makeop();
    linop_adj = inn_linop.makeadj();
    linop_toiter = @(x) linop_adj(linop_op(x));
    eigsOpt.issym = 1;
    eigsOpt.tol = 1e-3;
    sqnorm_linop = eigs(linop_toiter, F.n, 1, 'LM', eigsOpt);
else
    sqnorm_linop = opt.sqOpNorm;
end

for it = 1:opt.maxit

    % define smooth term of the inner augmented Lagrangian subproblem
    hgamma = moreauEnvelope(h, 1/r);
    inn_f = separableSum({hgamma, f}, [F.m, F.n]);
    inn_aff = {inn_linop, [y/r; zeros(F.n, 1)]};
    if isfield(f, 'L'), inn_opt.Lf = f.L + sqnorm_linop*r; end

    % solve subproblem (warm start)
    inn_out = forbes(inn_f, g, x, inn_aff, [], inn_opt);
    x = inn_out.x;
    tot_inn_it = tot_inn_it + inn_out.iterations;
    tot_ops = OpsSum(tot_ops, inn_out.operations);

    % compute next dual iterate
    callhgamma = hgamma.makef();
    [~, y1] = callhgamma(callF(x) + y/r);
    res(1,it) = norm(y1-y)/r;
    
    % display info
    if opt.display == 1
        fprintf('.');
    elseif opt.display >= 2
        fprintf('%6d %7.4e %7.4e %7.4e %10d\n', it, res(1,it), r, inn_opt.tol, inn_out.iterations);
    end

    if res(1,it) <= eta
        if res(1,it) <= opt.tol && inn_out.residual(end)<=opt.tol
            break
        else
            y = y1;
            eta = eta/r^0.9;
            inn_opt.tol = inn_opt.tol/10;
        end
    else
        r = 10*r;
        eta = 1/r^0.1;
        inn_opt.tol = 1/r;
    end

    % stopping criterion
    if res(1, it) <= opt.tol
        break;
    end

end

if opt.display == 1
    fprintf('\n');
end

out.x = inn_out.x;
out.y = y;
out.iterations = it;
out.inner_iterations = tot_inn_it;
out.operations = tot_ops;
    
function opt = default_opt(opt)

if ~isfield(opt, 'display'), opt.display = 1; end
if ~isfield(opt, 'maxit'), opt.maxit = 100; end
if ~isfield(opt, 'tol'), opt.tol = 1e-6; end

function opt = default_inner_opt(opt)

if ~isfield(opt, 'display'), opt.display = 0; end
if ~isfield(opt, 'tol'), opt.tol = 1e-6; end
if ~isfield(opt, 'method'), opt.method = 'lbfgs-fpr'; end
% make sure the Lipschitz constant is not set (it cannot be known)
if isfield(opt, 'Lf'), opt = rmfield(opt, 'Lf'); end
