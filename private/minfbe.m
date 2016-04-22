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

function out = minfbe(prob, opt)

lsopt = ProcessLineSearchOptions(prob, opt);

% initialize output stuff
objective = zeros(1, opt.maxit);
ts = zeros(1, opt.maxit);
residual = zeros(1, opt.maxit);
msgTerm = '';

% initialize operations counter
ops = OpsInit();

% initialize stuff
x = prob.x0;
y = x;
gam = (1-opt.beta)/prob.Lf;

% initialize specific stuff for the methods
if opt.method == 2
    S = zeros(prob.n, opt.memory);
    Y = zeros(prob.n, opt.memory);
    YS = zeros(opt.memory, 1);
    LBFGS_col = 1;
    LBFGS_mem = 0;
end

% display header
if opt.display >= 2
    fprintf('%6s%11s%11s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', '||dir||', 'slope', 'tau');
end

if opt.toRecord
    record = [];
end

skipCount = 0;
flagChangedGamma = 0;

cache_current = CacheInit(prob, y, gam);

t0 = tic();

for it = 1:opt.maxit
    
    % compute FBE
    [cache_current, ops1] = CacheFBE(cache_current, gam);
    ops = OpsSum(ops, ops1);
    
%     % backtracking on gamma
%     cache_z = CacheInit(prob, cache_current.z, gam);
%     if prob.unknownLf || opt.adaptive
%         [cache_z, ops1] = CacheEvalf(cache_z);
%         ops = OpsSum(ops, ops1);
%         fz = cache_z.fx;
%         while fz + cache_current.gz + lsopt.beta/(2*gam)*cache_current.normdiff^2 > cache_current.FBE
%             prob.Lf = 2*prob.Lf; gam = gam/2;
%             flagChangedGamma = 1;
%             [cache_current, ops1] = CacheFBE(cache_current, gam);
%             ops = OpsSum(ops, ops1);
%             cache_z = CacheInit(prob, cache_current.z, gam);
%             [cache_z, ops1] = CacheEvalf(cache_z);
%             ops = OpsSum(ops, ops1);
%             fz = cache_z.fx;
%         end
%     end

    % store initial cache
    if it == 1
        cache_0 = cache_current;
    end

    % trace stuff
    ts(1, it) = toc(t0);
    residual(1, it) = norm(cache_current.diff, 'inf')/gam;
    objective(1, it) = cache_current.FBE;
    if opt.toRecord
        record = [record, opt.record(prob, it, gam, cache_0, cache_current, ops)];
    end

    % check for termination
    if isnan(cache_current.normdiff)
        msgTerm = 'something went wrong';
        flagTerm = 1;
        break;
    end
    if residual(1, it) <= 100*eps % cache_current.normdiff/(1+norm(cache_current.x)) <= 100*eps
        msgTerm = 'reached optimum (fpr close to eps)';
        flagTerm = 0;
        break;
    end
    if ~flagChangedGamma
        if ~opt.customTerm
            if residual(1, it) <= opt.tol
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

    % compute search direction and slope
    switch opt.method
        case {1, 6} % STEEPEST DESCENT and BARZILAI-BORWEIN
            dir = -cache_current.gradFBE;
        case 2 % L-BFGS
            if it == 1 || flagChangedGamma
                dir = -cache_current.gradFBE; % use steepest descent direction initially
                LBFGS_col = 0; % last column of Sk, Yk that was filled in
                LBFGS_mem = 0; % current memory of the method
            else
                %%% x' - x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Sk = cache_current.x - cache_previous.x;
                Yk = cache_current.gradFBE - cache_previous.gradFBE;
                YSk = Yk'*Sk;
                if YSk > 0
                    LBFGS_col = 1+mod(LBFGS_col, opt.memory);
                    LBFGS_mem = min(LBFGS_mem+1, opt.memory);
                    S(:,LBFGS_col) = Sk;
                    Y(:,LBFGS_col) = Yk;
                    YS(LBFGS_col) = YSk;
                else
                    skipCount = skipCount+1;
                end
                if LBFGS_mem > 0
                    H = YS(LBFGS_col)/(Y(:,LBFGS_col)'*Y(:,LBFGS_col));
                    dir = LBFGS(S, Y, YS, H, -cache_current.gradFBE, int32(LBFGS_col), int32(LBFGS_mem));
                else
                    dir = -cache_current.gradFBE;
                end
            end
        case 7 % BFGS
            if it == 1 || flagChangedGamma
                dir = -cache_current.gradFBE;
                R = eye(prob.n);
            else
                Sk = cache_current.x - cache_previous.x;
                Yk = cache_current.gradFBE - cache_previous.gradFBE;
                YSk = Yk'*Sk;
                Bs = R'*(R*Sk);
                sBs = Sk'*Bs;
                if YSk > 0
                    R = cholupdate(cholupdate(R,Yk/sqrt(YSk)),Bs/sqrt(sBs),'-');
                else
                    skipCount = skipCount+1;
                end
                dir = -linsolve(R,linsolve(R,cache_current.gradFBE,opt.optsL),opt.optsU);
            end
        case 3 % CG-DESCENT
            if it == 1 || flagChangedGamma
                dir = -cache_current.gradFBE; % Initially use steepest descent direction
            else
                yy = cache_current.gradFBE-cache_previous.gradFBE;
                dy = dir'*yy;
                lambda = 1; %Hager-Zhang proposed lambda = 2 but Kou, Dai found that lambda = 1 is more efficient
                %                 lambda = 2-(dir'*yy)/((dir'*dir)*(yy'*yy));
                beta = ((yy-lambda*dir*(yy'*yy)/dy)'*cache_current.gradFBE)/dy;
                etak = -1/(norm(dir)*min(0.01,norm(cache_current.gradFBE)));
                beta = max(beta,etak);
                dir = -cache_current.gradFBE + beta*dir;
                if dir'*cache_current.gradFBE >= 0 % restart if not descent direction
                    dir = -cache_current.gradFBE;
                    skipCount = skipCount+1;
                end
            end
        case 4 % CG-PRP
            if it == 1 || flagChangedGamma
                dir = -cache_current.gradFBE; % Initially use steepest descent direction
            else
                yy = cache_current.gradFBE - cache_previous.gradFBE;
                beta = max((cache_current.gradFBE'*yy)/(cache_previous.gradFBE'*cache_previous.gradFBE),0);
                dir = -cache_current.gradFBE + beta*dir;
                if dir'*cache_current.gradFBE >= 0 % restart if not descent direction
                    dir = -cache_current.gradFBE;
                    skipCount = skipCount+1;
                end
            end
        case 5 % CG-DYHS
            if it == 1 || flagChangedGamma
                dir = -cache_current.gradFBE; % Initially use steepest descent direction
            else
                yy = cache_current.gradFBE - cache_previous.gradFBE;
                betaDY = (cache_current.gradFBE'*cache_current.gradFBE)/(dir'*yy);
                betaHS = (cache_current.gradFBE'*yy)/(dir'*yy);
                beta = max(0,min(betaHS,betaDY));
                dir = -cache_current.gradFBE + beta*dir;
                if dir'*cache_current.gradFBE >= 0 % restart if not descent direction
                    dir = -cache_current.gradFBE;
                    skipCount = skipCount+1;
                end
            end
        otherwise
            error('search direction not implemented')
    end
    
    slope = cache_current.gradFBE'*dir;

    % precompute stuff for the line search
    [cache_current, ops1] = CacheLineSearch(cache_current, dir);
    ops = OpsSum(ops, ops1);

    % set initial guess for the step length
    switch opt.method
        case {2, 7} % (L-)BFGS
            tau0 = 1.0;
        case 6 % BARZILAI-BORWEIN
            if it == 1 || flagChangedGamma
                tau0 = 1.0/norm(cache_current.gradFBE, inf);
            else
                Sk = cache_current.x-cache_previous.x;
                Yk = cache_current.gradFBE-cache_previous.gradFBE;
                tau0 = (Sk'*Sk)/(Sk'*Yk);
            end
        otherwise
            if it == 1 || flagChangedGamma
                xinf = norm(cache_current.x,inf);
                if xinf ~= 0
                    tau0 = lsopt.psi0*xinf/norm(cache_current.gradFBE, inf); % g is the gradient at x
                elseif cache_current.FBE ~= 0
                    % Check changes in Version 5.2 (3). psi0 is set equal
                    % to 0 when x=0;
                    tau0 = lsopt.psi0*abs(cache_current.FBE)/(cache_current.gradFBE'*cache_current.gradFBE);
                else
                    tau0 = 1;
                end
            else
                %                             lsopt.tau0 = lsopt.psi2*tau;
                %                             lsopt.tau0 = 2*(cache.FBE - FBE_old)/slope;
                tau0 = -2*max(cache_previous.FBE-cache_current.FBE, 10*eps)/slope;% Fletcher, pp. 38
                if lsopt.quadStep
                    tp = tau*lsopt.psi1;
                    [cache_tau, ops1] = DirFBE(cache_current, tp, 1);
                    ops = OpsSum(ops, ops1);
                    if cache_tau.FBE <= cache_current.FBE
                        % quadratic interpolation
                        q = cache_tau.FBE - cache_current.FBE - tp*slope;
                        if q > 0
                            tau0 = -(slope*tp^2)/(2*q);
                        end
                    end
                end
            end
    end

    % perform line search
    switch opt.linesearch
        case 0 % NO LINE SEARCH - UNIT STEPSIZE
            [cache_tau, ops1] = CacheFBE(prob, gam, cache_current.x+dir);
            tau = 1.0;
            info = 0;
        case 1 % ARMIJO LINE SEARCH WITH CUBIC INTERPOLATION
            [tau, cache_tau, cache_tau1, ops1, info] = ArmijoLS(cache_current, slope, tau0, lsopt);
        case 3 % LEMARECHAL (WOLF CONDITIONS)
            [tau, cache_tau, cache_tau1, ops1, info] = LemarechalLS(cache_current, slope, tau0, lsopt);
        case 7 % SIMPLE BACKTRACKING (NON-INCREASING LINE-SEARCH)
            [tau, cache_tau, cache_tau1, ops1, info] = BacktrackingLS(cache_current, tau0, lsopt);
        otherwise
            error('line search not implemented')
    end
    ops = OpsSum(ops, ops1);

    % check for line search fails
    if info > 0 && ~opt.global && ~opt.fast
        flagTerm = 2;
        msgTerm = strcat('line search failed at iteration', num2str(it));
    end

    % prepare next iteration, store current solution
    flagChangedGamma = 0;
    if info == -1 % gam was too large
        cache_previous = cache_current;
        prob.Lf = prob.Lf*2; gam = gam/2;
        flagChangedGamma = 1;
    elseif info > 0 % line-search failed
        cache_previous = cache_current;
        cache_current = CacheInit(prob, cache_current.z, gam);
    elseif opt.global % globalized line-search method
        cache_previous = cache_current;
        if ~isempty(cache_tau1)
            cache_current = cache_tau1;
        else
            cache_current = CacheInit(prob, cache_tau.z, gam);
        end
    else % basic line-search method
        cache_previous = cache_current;
        cache_current = cache_tau;
    end

    % display stuff
    if opt.display == 1
        PrintProgress(it);
    elseif opt.display >= 2
        fprintf('%6d %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %d\n', it, gam, residual(1,it), objective(1,it), norm(dir), slope, tau, info);
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
if info ~= 0
    out.x = cache_current.z;
else
    out.x = cache_tau.z;
end
out.iterations = it;
out.operations = ops;
out.residual = residual(1, 1:it);
out.objective = objective(1, 1:it);
out.ts = ts(1, 1:it);
out.prob = prob;
out.opt = opt;
out.gam = gam;
if opt.toRecord
    out.record = record;
end
out.skip = skipCount;

