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

function out = minfbe2(prob, opt)

lsopt = ProcessLineSearchOptions(prob, opt);

% initialize output stuff
objective = zeros(1, opt.maxit);
ts = zeros(1, opt.maxit);
residual = zeros(1, opt.maxit);
record = [];

% initialize operations counter
ops = OpsInit();

% initialize stuff
x = prob.x0;
w = x;
gam = (1-opt.beta)/prob.Lf;

% initialize specific stuff for the methods
if opt.method == 2
    S = zeros(prob.n, opt.memory);
    Y = zeros(prob.n, opt.memory);
    YS = zeros(opt.memory, 1);
end

% display header
if opt.display >= 2
    fprintf('%6s%11s%11s%11s%11s%11s%11s\n', ...
        'iter', 'gamma', 'optim.', 'object.', '||dir||', 'slope', 'tau');
end

cntSkip = 0;
flagLS = 0;
flagChangedGamma = 1;
flagTerm = 1;
msgTerm = 'exceeded maximum iterations';

cache_x = CacheInit(prob, x, gam);
cache_w = cache_x;

t0 = tic();

for it = 1:opt.maxit
    
    [cache_w, ops1] = CacheFBE(cache_w, gam);
    ops = OpsSum(ops, ops1);
    
    % store initial cache
    if it == 1
        cache_0 = cache_w;
    end
    
    cache_Tw = CacheInit(prob, cache_w.z, gam);
    
    if prob.unknownLf || opt.adaptive
        [cache_Tw, ops1] = CacheEvalf(cache_Tw);
        ops = OpsSum(ops, ops1);
        phi_Tw = cache_Tw.fx + cache_w.gz;
        while phi_Tw + opt.beta/(2*gam)*cache_w.normdiff^2 > cache_w.FBE
            if flagChangedGamma, prob.Lf = 2*prob.Lf; gam = gam/2; end
            flagChangedGamma = 1;
            [cache_w, ops1] = CacheFBE(cache_x, gam); % w <-- x
            ops = OpsSum(ops, ops1);
            cache_Tw = CacheInit(prob, cache_w.z, gam);
            [cache_Tw, ops1] = CacheEvalf(cache_Tw);
            ops = OpsSum(ops, ops1);
            phi_Tw = cache_Tw.fx + cache_w.gz;
        end
    end

    % compute FBE and its gradient
    [cache_xnew, ops1] = CacheFBE(cache_Tw, gam); % xnew <-- Tw
    ops = OpsSum(ops, ops1);
    [cache_xnew, ops1] = CacheGradFBE(cache_xnew, gam);
    ops = OpsSum(ops, ops1);
    
    % trace stuff
    ts(it) = toc(t0);
    residual(it) = norm(cache_xnew.diff, 'inf')/gam;
    objective(it) = cache_xnew.FBE;
    if opt.toRecord
        record = [record, opt.record(prob, it, gam, cache_0, cache_xnew, ops)];
    end
    
    % stopping criterion
    if ~flagChangedGamma
        if ~opt.customTerm
            if StoppingCriterion(cache_xnew, opt.tol)
                msgTerm = 'reached optimum (up to tolerance)';
                flagTerm = 0;
                break;
            end
        else
            flagStop = opt.term(prob, it, gam, cache_0, cache_xnew, ops);
            if (prob.unknownLf == 0 || it > 1) && flagStop
                msgTerm = 'reached optimum (custom criterion)';
                flagTerm = 0;
                break;
            end
        end
    end

    % compute search direction and slope
    switch opt.method
        case 2 % BFGS
            if flagChangedGamma || flagLS
                direction = -cache_xnew.gradFBE;
                R = eye(prob.n);
            else
                Sk = cache_xnew.x - cache_x.x;
                Yk = cache_xnew.gradFBE - cache_x.gradFBE;
                YSk = Yk'*Sk;
                Bs = R'*(R*Sk);
                sBs = Sk'*Bs;
                if YSk > 0
                    R = cholupdate(cholupdate(R,Yk/sqrt(YSk)),Bs/sqrt(sBs),'-');
                else
                    cntSkip = cntSkip+1;
                end
                direction = -linsolve(R,linsolve(R,cache_xnew.gradFBE,opt.optsL),opt.optsU);
            end
        case 3 % L-BFGS
            if flagChangedGamma || flagLS
                direction = -cache_xnew.gradFBE; % use steepest descent direction initially
                LBFGS_col = 0; % last column of Sk, Yk that was filled in
                LBFGS_mem = 0; % current memory of the method
            else
                %%% x' - x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Sk = cache_xnew.x - cache_x.x;
                Yk = cache_xnew.gradFBE - cache_x.gradFBE;
                YSk = Yk'*Sk;
                if YSk > 0
                    LBFGS_col = 1+mod(LBFGS_col, opt.memory);
                    LBFGS_mem = min(LBFGS_mem+1, opt.memory);
                    S(:,LBFGS_col) = Sk;
                    Y(:,LBFGS_col) = Yk;
                    YS(LBFGS_col) = YSk;
                else
                    cntSkip = cntSkip+1;
                end
                if LBFGS_mem > 0
                    H = YS(LBFGS_col)/(Y(:,LBFGS_col)'*Y(:,LBFGS_col));
                    direction = LBFGS(S, Y, YS, H, -cache_xnew.gradFBE, int32(LBFGS_col), int32(LBFGS_mem));
                else
                    direction = -cache_xnew.gradFBE;
                end
            end
        otherwise
            error('search direction not implemented');
    end
    
    slope = cache_xnew.gradFBE'*direction;

    % set initial guess for the step length
    switch opt.method
        case {2, 3} % (L-)BFGS
            tau0 = 1.0;
        otherwise
            error('search direction not implemented');
    end
    
    % perform line search
    [cache_xnew, ops1] = CacheLineSearch(cache_xnew, direction);
    ops = OpsSum(ops, ops1);
    tau = tau0;
    flagLS = 1;
    f0 = cache_xnew.FBE;
    for i = 1:lsopt.nLS
        [cache_wnew, ops1] = DirFBE(cache_xnew, tau, 1);
        ops = OpsSum(ops, ops1);
        ft = cache_wnew.FBE;
        if ft <= f0
            flagLS = 0;
            break;
        end
        tau = 0.5*tau;
        if tau <= lsopt.progTol
            flagLS = 2;
            break
        end
    end

    % prepare next iteration
    flagChangedGamma = 0;
    cache_x = cache_xnew;
    if flagLS == 0, cache_w = cache_wnew;
    else cache_w = cache_xnew; end
    
    % display stuff
    if opt.display == 1
        PrintProgress(it);
    elseif opt.display >= 2
        fprintf('%6d %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %d\n', ...
            it, gam, residual(1,it), objective(1,it), norm(direction), slope, tau, flagLS);
    end

end

if opt.display == 1
    PrintProgress(it, flagTerm);
end

% pack up results
out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
out.x = cache_w.z;
out.iterations = it;
out.operations = ops;
out.residual = residual(1, 1:it);
out.objective = objective(1, 1:it);
out.record = record;
out.ts = ts(1, 1:it);
out.prob = prob;
out.opt = opt;
out.gam = gam;
out.skip = cntSkip;
