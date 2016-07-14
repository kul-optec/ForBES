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

function opt = ProcessOptions(opt)

% fill in missing options with defaults
if ~isfield(opt, 'tol') || isempty(opt.tol), opt.tol = 1e-8; end
if ~isfield(opt, 'term') || isempty(opt.term), opt.customTerm = false;
else opt.customTerm = true; end
if ~isfield(opt, 'record') || isempty(opt.record), opt.toRecord = false;
else opt.toRecord = true; end
if ~isfield(opt, 'maxit') || isempty(opt.maxit), opt.maxit = 10000; end
if ~isfield(opt, 'beta') || isempty(opt.beta), opt.beta = 0.05; end
if ~isfield(opt, 'solver') || isempty(opt.solver), opt.solver = 'zerofpr'; end
if ~isfield(opt, 'adaptive') || isempty(opt.adaptive), opt.adaptive = 0; end
if ~isfield(opt, 'display') || isempty(opt.display), opt.display = 2; end
if ~isfield(opt, 'useHessian') || isempty(opt.useHessian), opt.useHessian = 0; end

opt.name = opt.solver;

% translate labels into integer codes
switch opt.solver
    case 'fbs'
        opt.solverID = 1;
    case 'minfbe'
        opt.solverID = 2;
    case 'zerofpr'
        opt.solverID = 3;
        if ~isfield(opt, 'qnopt') || isempty(opt.qnopt), opt.qnopt = 1; end
    case 'minfbe2'
        opt.solverID = 4;
    otherwise
        error('unknown solver');
end

solver2variant = {'fast', 'global', '', 'global'};
if ~isfield(opt, 'variant') || isempty(opt.variant)
    opt.variant = solver2variant{opt.solverID};
end

opt.name = strcat(opt.name, ',' , opt.variant);

opt.fast = 0;
opt.global = 0;
switch opt.variant
    case 'global'
        opt.global = 1;
    case 'fast'
        opt.fast = 1;
end

solver2method = {'none', 'lbfgs', 'lbfgs', 'lbfgs'};
if ~isfield(opt, 'method') || isempty(opt.method)
    opt.method = solver2method{opt.solverID};
end

opt.name = strcat(opt.name, ', ', opt.method);

switch opt.method
    case 'none'
        opt.methodID = 0;
    case 'sd'
        opt.methodID = 1;
    case 'bfgs'
        opt.methodID = 2;
        opt.optsL.UT = true;
        opt.optsL.TRANSA = true;
        opt.optsU.UT = true;
    case 'lbfgs'
        opt.methodID = 3;
        opt.memory = 10;
    case 'cg-desc'
        opt.methodID = 4;
    case 'cg-prp'
        opt.methodID = 5;
    case 'cg-dyhs'
        opt.methodID = 6;
    case 'bb'
        opt.methodID = 7;
    case 'broyden'
        opt.methodID = 8;
        if ~isfield(opt, 'bopt') || isempty(opt.bopt), opt.bopt = 0; end
    case 'lbroyden'
        opt.methodID = 9;
        opt.memory = 5;
        if ~isfield(opt, 'bopt') || isempty(opt.bopt), opt.bopt = 0; end
    case 'rbroyden'
        opt.methodID = 10;
        opt.memory = 5;
        if ~isfield(opt, 'bopt') || isempty(opt.bopt), opt.bopt = 0; end
    otherwise
        error('unknown method');
end

% the default line-search depends on solver, method and variant
if opt.methodID > 0 && (~isfield(opt, 'linesearch') || isempty(opt.linesearch))
    if (opt.solverID == 2 || opt.solverID == 4) && opt.global == 0 % minfbe/minfbe2, classical line search method
        method2linesearch = {'armijo', 'lemarechal', 'lemarechal', 'lemarechal', 'lemarechal', 'lemarechal', 'nonmonotone-armijo', 'lemarechal', 'lemarechal'};
        opt.linesearch = method2linesearch{opt.methodID};
    elseif (opt.solverID == 2 || opt.solverID == 4) && opt.global == 1 % minfbe/minfbe2, classical line search method
        method2linesearch = repmat({'backtracking'}, 1, 9);
        opt.linesearch = method2linesearch{opt.methodID};
    elseif (opt.solverID == 3)
        method2linesearch = repmat({'backtracking'}, 1, 9);
        opt.linesearch = method2linesearch{opt.methodID};
    else
        opt.linesearch = 'none';
    end
elseif opt.methodID == 0
    opt.linesearch = 'none';
end

opt.name = strcat(opt.name, ', ', opt.linesearch);

switch opt.linesearch
    case 'none'
        opt.linesearchID = 0;
    case 'backtracking'
        opt.linesearchID = 1;
    case 'backtracking-nm'
        opt.linesearchID = 2;
    case 'backtracking-armijo'
        opt.linesearchID = 3;
    case 'lemarechal'
        opt.linesearchID = 4;
    case 'hager-zhang'
        opt.linesearchID = 5;
    case 'more-thuente'
        opt.linesearchID = 6;
    case 'fletcher'
        opt.linesearchID = 7;
    otherwise
        error('unknown line search');
end

opt.processed = true;
