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

function out = zerofpr(prob, opt)

% initialize line-search options

lsopt = ProcessLineSearchOptions(prob, opt);
lsopt.testGamma = 0;

% initialize operations counter

ops = OpsInit();

% initialize gamma and sigma

gam = (1-opt.beta)/prob.Lf;
sig = (1-gam*prob.Lf)/(4*gam);

% display header

if opt.display >= 2
    fprintf('%6s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', 'tau');
end

cacheDir.cntSkip = 0;

alpha = 0.5;
flag = -1; % to track what happened at every iteration
flagTerm = 0;

MAXIMUM_Lf = 1e15;
MINIMUM_tau = 1e-15;
MINIMUM_d = 1e-15;

t0 = tic();

cache_x = CacheInit(prob, prob.x0, gam);
[cache_x, ops1] = CacheProxGradStep(cache_x, gam);
ops = OpsSum(ops, ops1);
cache_0 = cache_x;

for it = 1:opt.maxit
    
    hasGammaChanged = 0;
    
    % backtracking on gamma
    
    if prob.unknownLf || opt.adaptive
        [isGammaOK, cache_x, cache_xbar, ops1] = CheckGamma(cache_x, gam, opt.beta);
        ops = OpsSum(ops, ops1);
        while ~isGammaOK
            prob.Lf = 2*prob.Lf; gam = gam/2; sig = 2*sig;
            hasGammaChanged = 1;
            [isGammaOK, cache_x, cache_xbar, ops1] = CheckGamma(cache_x, gam, opt.beta);
            ops = OpsSum(ops, ops1);
        end
    else
        cache_xbar = CacheInit(prob, cache_x.z, gam);
    end
    
    if prob.Lf >= MAXIMUM_Lf
        msgTerm = ['estimate for Lf became too large: ', num2str(prob.Lf)];
        flagTerm = 1;
        break;
    end
    
    % trace stuff
    
    ts(1, it) = toc(t0);
    residual(1, it) = norm(cache_x.FPR, 'inf')/gam;
    if opt.toRecord
        record(:, it) = opt.record(prob, it, gam, cache_0, cache_x, ops);
    end
    
    % compute FBE at current point
    % this should count zero operations if gamma hasn't changed
    
    [cache_x, ops1] = CacheFBE(cache_x, gam);
    ops = OpsSum(ops, ops1);
    
    objective(1,it) = cache_x.FBE;
    
    % check for termination
    
    if isnan(cache_x.normFPR)
        msgTerm = 'something went wrong';
        flagTerm = 1;
        break;
    end
    if ~hasGammaChanged
        if ~opt.customTerm 
            if StoppingCriterion(cache_x, opt.tol)
                msgTerm = 'reached optimum (up to tolerance)';
                flagTerm = 0;
                break;
            end
        else
            flagStop = opt.term(prob, it, gam, cache_0, cache_x, ops);
            if (prob.unknownLf == 0 || it > 1) && flagStop
                msgTerm = 'reached optimum (custom criterion)';
                flagTerm = 0;
                break;
            end
        end
    end
    
    % select a direction
    
    [cache_xbar, ops1] = CacheProxGradStep(cache_xbar, gam);
    ops = OpsSum(ops, ops1);
    
    % compute pair (s, y) for quasi-Newton updates
    
    if it > 1 && ~hasGammaChanged
        sk = cache_x.x - cache_previous.x;
        yk = cache_x.FPR - cache_previous.FPR;
%         fprintf('%7.4e %7.4e %7.4e\n', norm(yk - opt.G_star(gam, sk))/norm(sk), norm(yk - opt.H_star(gam, sk))/norm(sk), norm(yk - opt.G_sym_star(gam, sk))/norm(sk));
    else
        sk = [];
        yk = [];
    end

    % compute search direction and slope
    
    [dir, cacheDir] = ComputeDir(prob, opt, it, hasGammaChanged, sk, yk, ...
        cache_xbar.FPR, cacheDir);
    cacheDir.prev_v = cache_xbar.FPR;
    
    % set initial guess for the step length
    
    tau0 = ComputeTau0(prob, opt, sk, yk, dir);
    
    % perform line search
    
    switch opt.linesearchID
        case 1 % backtracking
            ref = cache_x.FBE - sig*cache_x.normFPR^2;
            [tau, cache_tau, ~, ops1, ~] = BacktrackingLS(cache_xbar, dir, tau0, lsopt, ref);
        case 2 % backtracking (nonmonotone)
            if it == 1 || hasGammaChanged
                Q = 1;
                C = cache_x.FBE;
            else
                newQ = lsopt.eta*Q+1;
                C = (lsopt.eta*Q*C + cache_x.FBE)/newQ;
                Q = newQ;
            end
            ref = C - sig*cache_x.normFPR^2;
            [tau, cache_tau, ~, ops1, ~] = BacktrackingLS(cache_xbar, dir, tau0, lsopt, ref);
        otherwise
            error('line search not implemented')
    end
    cacheDir.prev_tau = tau;
    
    % update iterates
    
    if opt.qnopt == 1
        cache_previous = cache_xbar; % s = x_{k+1}-\bar{x}_{k}, y = r_{k+1}-\bar{r}_{k}
    elseif opt.qnopt == 2
        cache_previous = cache_x; % s = x_{k+1}-x_{k}, y = r_{k+1}-r_{k}
    end
    
    cache_x = cache_tau;
    ops = OpsSum(ops, ops1);

    if flagTerm == 1
        break;
    end

    % display stuff
    
    if opt.display == 1
        PrintProgress(it);
    elseif opt.display >= 2
        fprintf('%6d %7.4e %7.4e %7.4e %7.4e\n', it, gam, residual(1,it), objective(1,it), tau);
    end

end

if it == opt.maxit
    msgTerm = 'exceeded maximum iterations';
    flagTerm = 1;
end

if opt.display == 1
    PrintProgress(it, flagTerm);
end

% pack up results
out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
out.x = cache_x.z;
out.iterations = it;
out.operations = ops;
out.residual = residual(1, 1:it);
out.objective = objective(1, 1:it);
out.ts = ts(1, 1:it);
if opt.toRecord, out.record = record; end
out.prob = prob;
out.opt = opt;
out.gam = gam;
