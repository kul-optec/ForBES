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

function out = minfbe(prob, opt)

% initialize line-search options

lsopt = ProcessLineSearchOptions(prob, opt);

% initialize operations counter

ops = OpsInit();

% initialize gamma

gam = (1-opt.beta)/prob.Lf;

% display header

if opt.display >= 2
    fprintf('%6s%11s%11s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', '||dir||', 'slope', 'tau');
end

cacheDir.cntSkip = 0;
hasGammaChanged = 0;

cache_current = CacheInit(prob, prob.x0, gam);

t0 = tic();

for it = 1:opt.maxit

    % compute FBE

    [cache_current, ops1] = CacheFBE(cache_current, gam);
    ops = OpsSum(ops, ops1);

    % store initial cache

    if it == 1
        cache_0 = cache_current;
    end

    % trace stuff

    ts(1, it) = toc(t0);
    residual(1, it) = norm(cache_current.FPR, 'inf')/gam;
    objective(1, it) = cache_current.FBE;
    if opt.toRecord
        record(:, it) = opt.record(prob, it, gam, cache_0, cache_current, ops);
    end

    solution = cache_current.z;

    % check for termination

    if isnan(cache_current.normFPR)
        msgTerm = 'something went wrong';
        flagTerm = 1;
        break;
    end
    if ~hasGammaChanged
        if ~opt.customTerm
            if StoppingCriterion(cache_current, opt.tol)
                msgTerm = 'reached optimum (up to tolerance)';
                flagTerm = 0;
                break;
            end
        else
            flagStop = opt.term(prob, it, gam, cache_0, cache_current, ops);
            if (prob.unknownLf == 0 || it > 1) && flagStop
                msgTerm = 'reached optimum (custom criterion)';
                flagTerm = 0;
                break;
            end
        end
    end

    % compute gradient of the FBE

    [cache_current, ops1] = CacheGradFBE(cache_current, gam);
    ops = OpsSum(ops, ops1);

    % compute pair (s, y) for quasi-Newton updates

    if it > 1
        sk = cache_current.x - cache_previous.x;
        yk = cache_current.gradFBE - cache_previous.gradFBE;
    else
        sk = [];
        yk = [];
    end

    % compute search direction and slope

    [dir, cacheDir] = ComputeDir(prob, opt, it, hasGammaChanged, sk, yk, ...
        cache_current.gradFBE, cacheDir);
    slope = cache_current.gradFBE'*dir;

    % set initial guess for the step length

    tau0 = ComputeTau0(prob, opt, sk, yk, dir);

    % perform line search

    switch opt.linesearchID
        case 1 % backtracking
            [tau, cache_tau, cache_tau1, ops1, flagLS] = BacktrackingLS(cache_current, dir, tau0, lsopt);
        case 2 % backtracking (nonmonotone)
            if it == 1 || hasGammaChanged
                Q = 1;
                C = cache_current.FBE;
            else
                newQ = lsopt.eta*Q+1;
                C = (lsopt.eta*Q*C + cache_current.FBE)/newQ;
                Q = newQ;
            end
            [tau, cache_tau, cache_tau1, ops1, flagLS] = BacktrackingLS(cache_current, dir, tau0, lsopt, C);
        case 3 % backtracking (Armijo condition)
            ref = cache_current.FBE; % f(0)
            lin = lsopt.delta*slope; % delta f'(0)
            [tau, cache_tau, cache_tau1, ops1, flagLS] = BacktrackingLS(cache_current, dir, tau0, lsopt, ref, lin);
        case 4 % Lemarechal (Wolfe conditions)
            [tau, cache_tau, cache_tau1, ops1, flagLS] = LemarechalLS(cache_current, dir, slope, tau0, lsopt);
        otherwise
            error('line search not implemented')
    end
    ops = OpsSum(ops, ops1);

    % check for line search fails

    if flagLS > 0 && ~opt.global
        flagTerm = 2;
        msgTerm = strcat('line search failed at iteration', num2str(it));
        break;
    end

    % prepare next iteration, store current solution

    hasGammaChanged = 0;
    if flagLS == -1 % gam was too large
        cache_previous = cache_current;
        prob.Lf = prob.Lf*2; gam = gam/2;
        hasGammaChanged = 1;
        solution = cache_current.z;
    elseif flagLS > 0 % line-search failed
        cache_previous = cache_current;
        cache_current = CacheInit(prob, cache_current.z, gam);
        solution = cache_current.z;
    elseif opt.global % globalized line-search method
        cache_previous = cache_current;
        if ~isempty(cache_tau1)
            cache_current = cache_tau1;
            solution = cache_previous.z;
        else
            cache_current = CacheInit(prob, cache_tau.z, gam);
            solution = cache_tau.z;
        end
    else % basic line-search method
        cache_previous = cache_current;
        cache_current = cache_tau;
        solution = cache_tau.z;
    end

    % display stuff

    if opt.display == 1
        PrintProgress(it);
    elseif opt.display >= 2
        fprintf('%6d %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %d\n', it, gam, residual(1,it), objective(1,it), norm(dir), slope, tau, flagLS);
    end

end

if it == opt.maxit
    flagTerm = 1;
    msgTerm = 'exceeded maximum iterations';
end

if opt.display == 1
    PrintProgress(it, flagTerm);
end

% pack up results
out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
out.x = solution;
out.iterations = it;
out.operations = ops;
out.residual = residual(1, 1:it);
out.objective = objective(1, 1:it);
out.ts = ts(1, 1:it);
if opt.toRecord, out.record = record; end
out.prob = prob;
out.opt = opt;
out.gam = gam;
out.skip = cacheDir.cntSkip;
