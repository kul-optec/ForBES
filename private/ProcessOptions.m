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
    if ~isfield(opt, 'tol'), opt.tol = 1e-8; end
    if ~isfield(opt, 'term'), opt.customTerm = false;
    else opt.customTerm = true; end
    if ~isfield(opt, 'record'), opt.toRecord = false;
    else opt.toRecord = true; end
    if ~isfield(opt, 'maxit'), opt.maxit = 10000; end
    if ~isfield(opt, 'method'), opt.method = 'lbfgs'; end
    
    if ~isfield(opt, 'linesearch')
        switch opt.method
            case 'sd'
                opt.linesearch = 'armijo';
            case 'lbfgs'
                opt.linesearch = 'lemarechal';
            case 'bfgs'
                opt.linesearch = 'lemarechal';
            case 'cg-desc'
                opt.linesearch = 'lemarechal';
            case 'cg-prp'
                opt.linesearch = 'lemarechal';
            case 'cg-dyhs'
                opt.linesearch = 'lemarechal';
            case 'bb'
                opt.linesearch = 'nonmonotone-armijo';
            case 'fbs'
                opt.linesearch = 'none';
            case 'broyden-fpr'
                opt.linesearch = 'none';
            case 'bfgs-fpr'
                opt.linesearch = 'none';
            case 'lbfgs-fpr'
                opt.linesearch = 'none';
            otherwise
                error('unknown method');
        end
    end
    
    opt.name = strcat(opt.method, ', ', opt.linesearch);
    
    if ~isfield(opt, 'recache'), opt.recache = 100; end
    if ~isfield(opt, 'memory'), opt.memory = 10; end
    if ~isfield(opt, 'adaptive'), opt.adaptive = 0; end
    if ~isfield(opt, 'display'), opt.display = 2; end
    if ~isfield(opt, 'useHessian'), opt.useHessian = 0; end
    
    % translate labels into integer codes
    switch opt.method
        case 'sd'
            opt.method = 1;
        case 'lbfgs'
            opt.method = 2;
        case 'cg-desc'
            opt.method = 3;
        case 'cg-prp'
            opt.method = 4;
        case 'cg-dyhs'
            opt.method = 5;
        case 'bb'
            opt.method = 6;
        case 'bfgs'
            opt.method = 7;
            opt.optsL.UT = true; opt.optsL.TRANSA = true;
            opt.optsU.UT = true;
        case 'fbs'
            opt.method = 0;
        case 'broyden-fpr'
            opt.method = 11;
        case 'bfgs-fpr'
            opt.method = 12;
        case 'lbfgs-fpr'
            opt.method = 13;
        otherwise
            error('unknown method');
    end
    
    switch opt.linesearch
        case 'armijo'
            opt.linesearch = 1;
        case 'nonmonotone-armijo'
            opt.linesearch = 2;
        case 'lemarechal'
            opt.linesearch = 3;
        case 'hager-zhang'
            opt.linesearch = 4;
        case 'more-thuente'
            opt.linesearch = 5;
        case 'fletcher'
            opt.linesearch = 6;
        case 'none'
            opt.linesearch = 0;
        otherwise
            error('unknown line search');
    end
    
    if ~isfield(opt, 'variant')
        switch opt.method
            case 0 % FBS
                opt.variant = 'fast';
            otherwise
                opt.variant = 'global';
        end
    end
    
    switch opt.variant
        case 'basic'
            opt.fast = 0;
            opt.global = 0;
            opt.monotone = 0;
        case 'global'
            opt.fast = 0;
            opt.global = 1;
            opt.monotone = 0;
        case 'fast'
            opt.fast = 1;
            opt.global = 0;
            opt.monotone = 1;
        otherwise
            error('unknown variant');
    end
    
    opt.processed = true;
end
