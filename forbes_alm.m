%FORBES_ALM Solver for nonsmooth optimization problems.
%
%   FORBES_ALM(f, g, h, F, init, opt, inn_opt) solves problems of the form
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

function out = forbes_alm(f, g, h, F, init, opt, inn_opt)

if nargin < 6, opt = []; end
if nargin < 7, inn_opt = []; end

% fill-in missing options with defaults
opt = default_opt(opt);
inn_opt = default_inner_opt(inn_opt);

if opt.display >= 2
    fprintf('%6s%11s%11s%11s\n', 'iter', 'res', 'penalty', 'inner it');
end

y = init;
x = zeros(F.n,1);
r = 1; % this must be done properly
res = zeros(1,opt.maxit);
callF = F.makeop();

tot_inn_it = 0;
tot_ops = OpsInit();

max_tol = inn_opt.tol;

for it = 1:opt.maxit

    % define augmented Lagrangian subproblem
    hgamma = moreauEnvelope(h, 1/r);
    inn_f = separableSum({hgamma, f}, [F.m, F.n]);
    inn_aff = {stackOp({F, diagOp(F.n)}), [y/r; zeros(F.n, 1)]};

    % solve subproblem (warm start)
    inn_out = forbes(inn_f, g, x, inn_aff, [], inn_opt);
    x = inn_out.x;
    tot_inn_it = tot_inn_it + inn_out.iterations;
    tot_ops = OpsSum(tot_ops, inn_out.operations);

    % compute next dual iterate
    callhgamma = hgamma.makef();
    [~, y1] = callhgamma(callF(x) + y/r);
    res(1,it) = norm(y1-y)/r;
    inn_opt.tol = min(max_tol, 1e-2*res(1,it));

    % update dual variable
    y = y1;

    % display info
    if opt.display == 1
        fprintf('.');
    elseif opt.display >= 2
        fprintf('%6d %7.4e %7.4e %10d\n', it, res(1,it), r, inn_out.iterations);
    end

    % stopping criterion
    if res(1, it) <= opt.tol
        break;
    end

    % adjust penalty parameter r
    if it > 1 && res(1,it) >= 0.5*res(1,it-1)
        r = min(1e6, 2*r);
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
