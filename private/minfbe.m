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
    
    %% initialize output stuff
    objective = zeros(1, opt.maxit);
    ts = zeros(1, opt.maxit);
    residual = zeros(1, opt.maxit);
    msgTerm = '';
    
    %% initialize counters
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];

    %% initialize stuff
    x = prob.x0;
    y = x;
    recache = true;
    gam = SelectGamma(prob, opt);
    
    %% initialize specific stuff for the methods
    if opt.method == 2
        S = zeros(prob.n, opt.memory);
        Y = zeros(prob.n, opt.memory);
        YS = zeros(opt.memory, 1);
        LBFGS_col = 1;
        LBFGS_mem = 0;
    end
    if opt.linesearch == 2 % nonmonotone Armijo
        FBErefs = -inf*ones(lsopt.M, 1);
    end
    Q = 0; C = 0; % parameters for Hager-Zhang line search
    
    %% display stuff
    if opt.display >= 2
        fprintf('%6s%11s%11s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', '||dir||', 'slope', 'tau');
    end
    
    if opt.toRecord
        record = [];
    end
    
    skipCount = 0;
    rejCount = 0;
    flagChangedGamma = 0;
    
    t0 = tic();
    
    %% main iteration
    for it = 1:opt.maxit
        %% compute proximal gradient step and FBE (if needed)
        if recache || mod(it, opt.recache) == 0
            [cache_current, cnt1] = CacheFBE(prob, gam, y);
            cnt = cnt+cnt1;
            z = cache_current.z;
            % enforce monotonicity in the fast-monotone case
            if opt.fast && it > 1 && opt.monotone > 0 && ~flagChangedGamma
                if cache_current.FBE > cache_tau.FBE
                    y = cache_tau.z;
                    [cache_current, cnt1] = CacheFBE(prob, gam, y);
                    cnt = cnt+cnt1;
                    rejCount = rejCount+1;
                end
            end
        end
        
        %% store initial cache
        if it == 1
            cache_0 = cache_current;
        end
        
        %% trace stuff
        ts(1, it) = toc(t0);
        residual(1, it) = norm(cache_current.diff, inf)/gam;
        objective(1, it) = cache_current.FBE;
        if opt.toRecord
            record = [record, opt.record(prob, it, gam, cache_0, cache_current, cnt)];
        end
        
        %% check for termination
        if ~flagChangedGamma
            if ~opt.customTerm
                % From sec. 8.2.3.2 of Gill, Murray, Wright (1982).
                absFBE = abs(cache_current.FBE);
                flagStop = residual(1, it) <= nthroot(opt.tol, 3)*(1+absFBE) && ...
                        norm(cache_previous.x-cache_current.x, inf) < sqrt(opt.tol)*(1+norm(cache_current.x, inf)) && ...
                        cache_previous.FBE-cache_current.FBE < opt.tol*(1+absFBE);
                if residual(1, it) <= 10*sqrt(eps) || ((prob.unknownLf == 0 || it > 1) && flagStop)
                    msgTerm = 'reached optimum (up to tolerance)';
                    flagTerm = 0;
                    break;
                end
            else
                flagStop = opt.term(prob, it, gam, cache_0, cache_current, cnt);
                if (prob.unknownLf == 0 || it > 1) && flagStop
                    msgTerm = 'reached optimum (custom criterion)';
                    flagTerm = 0;
                    break;
                end
            end
        end
        
        %% compute gradient of the FBE
        [cache_current, cnt1] = CacheGradFBE(prob, gam, cache_current);
        cnt = cnt+cnt1;
        
        %% compute search direction and slope
        switch opt.method
            case {1, 6} % STEEPEST DESCENT and BARZILAI-BORWEIN
                dir = -cache_current.gradFBE;
            case 2 % L-BFGS
                if it == 1 || flagChangedGamma
                    dir = -cache_current.gradFBE; % use steepest descent direction initially
                    LBFGS_col = 1;
                    LBFGS_mem = 0;
                else
                    %%% x' - x
%                    Sk = cache_current.x - cache_previous.x;
%                    Yk = cache_current.gradFBE - cache_previous.gradFBE;
                    %%% other options (is this additional gradient eval needed?)
                    [cache_tau, cnt1] = CacheGradFBE(prob, gam, cache_tau);
                    cnt = cnt+cnt1;
                    %%% w - x
                    Sk = cache_tau.x - cache_previous.x;
                    Yk = cache_tau.gradFBE - cache_previous.gradFBE;
                    %%% x' - w
%                    Sk = cache_current.x - cache_tau.x;
%                    Yk = cache_current.gradFBE - cache_tau.gradFBE;
                    YSk = Yk'*Sk;
                    if YSk > 0;
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
%                elseif it == 2
%                    % Compute difference vectors
%                    Sk = cache_current.x - cache_previous.x;
%                    Yk = cache_current.gradFBE - cache_previous.gradFBE;
%                    YSk = Yk'*Sk;
%                    R = sqrt((Yk'*Yk)/(YSk))*eye(prob.n);
%                    %                     dir = -R\(R'\cache_current.gradFBE);
%                    dir = -cache_current.gradFBE./(diag(R).^2);
                else
                    %%% x' - x
%                     Sk = cache_current.x - cache_previous.x;
%                     Yk = cache_current.gradFBE - cache_previous.gradFBE;
                    %%% other options (is this additional gradient eval needed?)
                    [cache_tau, cnt1] = CacheGradFBE(prob, gam, cache_tau);
                    cnt = cnt+cnt1;
                    %%% w - x
                    Sk = cache_tau.x - cache_previous.x;
                    Yk = cache_tau.gradFBE - cache_previous.gradFBE;
                    %%% x' - w
%                     Sk = cache_current.x - cache_tau.x;
%                     Yk = cache_current.gradFBE - cache_tau.gradFBE;
                    YSk = Yk'*Sk;
                    Bs = R'*(R*Sk);
                    sBs = Sk'*Bs;
                    if YSk > 0
                        R = cholupdate(cholupdate(R,Yk/sqrt(YSk)),Bs/sqrt(sBs),'-');
                    else
                        skipCount = skipCount+1;
                    end
                    dir = -linsolve(R,linsolve(R,cache_current.gradFBE,opt.optsL),opt.optsU);
                    %                     dir = -R\(R'\cache_current.gradFBE);
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
        
        %% precompute stuff for the line search
        [cache_current, cnt1] = CacheLSData(prob, dir, cache_current);
        cnt = cnt+cnt1;

        %% set initial guess for the step length
        switch opt.method
            case 2 % L-BFGS
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
                    tau0 = -2*max(cache_previous.FBE-cache_current.FBE ,10*eps)/slope;% Fletcher, pp. 38
                    if lsopt.quadStep
                        tp = tau*lsopt.psi1;
                        [cache_tau, cnt1] = DirFBE(prob, gam, tp, cache_current, 1);
                        cnt = cnt+cnt1;
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

        %% perform the line search
        switch opt.linesearch
            case 0 % NO LINE SEARCH - UNIT STEPSIZE
                [cache_tau, cntLS] = CacheFBE(prob, gam, cache_current.x+dir);
                tau = 1.0;
                info = 0;
            case 1 % BACKTRACKING (ARMIJO COND.)
                [cache_tau, tau, cntLS, info] = ArmijoLS(prob, gam, cache_current, slope, tau0, lsopt);
            case 3 % LEMARECHAL (ARMIJO + WOLF COND.)
                [cache_tau, tau, cntLS, info] = LemarechalLS(prob, gam, cache_current, slope, tau0, lsopt);
            otherwise
                error('line search not implemented')
        end
        cnt = cnt+cntLS;

        %% check for line search fails
        if info > 0 && ~opt.global && ~opt.fast
            flagTerm = 2;
            msgTerm = strcat('line search failed at iteration', num2str(it));
        end
        
        flagChangedGamma = 0;
        if info == -1 % gam was too large
            y = cache_current.x;
            cache_previous = cache_current;
            recache = true;
            prob.Lf = prob.Lf*2;
            gam = SelectGamma(prob, opt);
            flagChangedGamma = 1;
        elseif info > 0
            y = cache_current.z;
            cache_previous = cache_current;
            recache = true;
        elseif opt.global
            y = cache_tau.z;
            cache_previous = cache_current;
            recache = true;
        elseif opt.fast
            theta1 = 2/(it+1);
            theta = 2/(it+2);
            v = x + (1/theta1)*(z-x);
            x = cache_tau.z;
            y = theta*v + (1-theta)*x;
            cache_previous = cache_current;
            recache = true;
        else
            y = cache_tau.x;
            cache_previous = cache_current;
            cache_current = cache_tau;
            recache = false;
        end
        
        %% display stuff
        if opt.display == 1
            if mod(it, 100) == 0
                fprintf('.');
            end
            if mod(it, 4000) == 0
                fprintf('\n');
            end
        elseif opt.display >= 2
            fprintf('%6d %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e %d\n', it, gam, residual(1,it), objective(1,it), norm(dir), slope, tau, info);
        end
        
    end
    
    if it == opt.maxit
        flagTerm = 1;
        msgTerm = 'exceeded maximum iterations';
    end
    
    %% pack up results
    out.name = opt.name;
    out.message = msgTerm;
    out.flag = flagTerm;
    out.x = cache_current.z;
    out.iterations = it;
    out.operations.cnt_f1 = cnt(1);
    out.operations.cnt_C1 = cnt(2);
    out.operations.cnt_f2 = cnt(4);
    out.operations.cnt_C2 = cnt(3);
    out.operations.cnt_g = cnt(5);
    out.operations.rej_extr = rejCount;
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
end

function gam = SelectGamma(prob, opt)
    if opt.method == 0 || opt.fast || opt.global
        gam = 1/prob.Lf;
    else
        gam = 0.95/prob.Lf;
    end
end

%% compute the FBE value and store reusable quantities
function [cache, cnt] = CacheFBE(prob, gam, x, cache)
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];
    f1x = 0; gradf1x = 0;
    f2x = 0; gradf2x = 0;
    cache.gam = gam;
    if nargin < 4
        cache.x = x;
        if prob.istheref1
            if prob.isthereC1
                if prob.isC1fun, C1x = prob.C1(cache.x);
                else C1x = prob.C1*cache.x; end
                cache.res1x = C1x - prob.d1;
                if prob.isQfun, cache.Qres1x = prob.Q(cache.res1x);
                else cache.Qres1x = prob.Q*cache.res1x; end
                if prob.isC1fun, gradf1x = prob.C1t(cache.Qres1x + prob.q);
                else gradf1x = prob.C1'*(cache.Qres1x + prob.q); end
                cnt(2) = cnt(2)+2;
            else
                cache.res1x = cache.x - prob.d1;
                if prob.isQfun, cache.Qres1x = prob.Q(cache.res1x);
                else cache.Qres1x = prob.Q*cache.res1x; end
                gradf1x = cache.Qres1x + prob.q;
            end
            cnt(1) = cnt(1)+1;
            cache.gradf1x = gradf1x;
            f1x = 0.5*(cache.res1x'*cache.Qres1x) + prob.q'*cache.res1x;
            cache.f1x = f1x;
        end
        if prob.istheref2
            if prob.isthereC2
                if prob.isC2fun, C2x = prob.C2(cache.x);
                else C2x = prob.C2*cache.x; end
                cache.res2x = C2x - prob.d2;
                if prob.useHessian
                    [f2x, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
                else
                    [f2x, gradf2res2x] = prob.callf2(cache.res2x);
                end
                if prob.isC2fun, gradf2x = prob.C2t(gradf2res2x);
                else gradf2x = prob.C2'*gradf2res2x; end
                cnt(3) = cnt(3)+2;
            else
                cache.res2x = cache.x - prob.d2;
                if prob.useHessian
                    [f2x, gradf2res2x, cache.Hessf2res2x] = prob.callf2(cache.res2x);
                else
                    [f2x, gradf2res2x] = prob.callf2(cache.res2x);
                end
                gradf2x = gradf2res2x;
            end
            cnt(4) = cnt(4)+1;
            cache.gradf2x = gradf2x;
        end
        if prob.istherelin
            cache.flinx = prob.l'*x;
            cache.fx = f1x + f2x + cache.flinx;
            cache.gradfx = gradf1x + gradf2x + prob.l;
        else
            cache.fx = f1x + f2x;
            cache.gradfx = gradf1x + gradf2x;
        end
    end
    y = cache.x - gam*cache.gradfx;
    [cache.z, cache.gz] = prob.callg(y, gam);
    cnt(5) = cnt(5)+1;
    cache.diff = cache.z-cache.x;
    sqnormdiff = cache.diff'*cache.diff;
    cache.normdiff = sqrt(sqnormdiff);
    cache.FBE = cache.fx + cache.gz + ...
        cache.gradfx'*cache.diff + ...
        (0.5/gam)*sqnormdiff;
end

%% compute the gradient of the FBE and store reusable quantities
function [cache, cnt] = CacheGradFBE(prob, gam, cache)
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];
    Hdiff = 0;
    if prob.istheref1
        if prob.isthereC1
            if prob.isC1fun, C1diff = prob.C1(cache.diff);
            else C1diff = prob.C1*cache.diff; end
            if prob.isQfun, QC1diff = prob.Q(C1diff);
            else QC1diff = prob.Q*C1diff; end
            if prob.isC1fun, C1tQC1diff = prob.C1t(QC1diff);
            else C1tQC1diff = prob.C1'*QC1diff; end
            cnt(2) = cnt(2)+2;
        else
            if prob.isQfun, C1tQC1diff = prob.Q(cache.diff);
            else C1tQC1diff = prob.Q*cache.diff; end
        end
        cnt(1) = cnt(1)+1;
        Hdiff = Hdiff + C1tQC1diff;
    end
    if prob.istheref2
        if prob.isthereC2
            if prob.isC2fun, C2diff = prob.C2(cache.diff);
            else C2diff = prob.C2*cache.diff; end
            cnt(3) = cnt(3)+1;
        else
            C2diff = cache.diff;
        end
        if prob.useHessian
            HC2diff = cache.Hessf2res2x*C2diff;
        else
            res2xepsdiff = cache.res2x + 1e-100i*C2diff;
            [~, gradf2res2xepsd] = prob.callf2(res2xepsdiff);
            cnt(4) = cnt(4)+1;
            HC2diff = imag(gradf2res2xepsd)/1e-100;
        end
        if prob.isthereC2
            if prob.isC2fun, Hdiff = Hdiff + prob.C2t(HC2diff);
            else Hdiff = Hdiff + (prob.C2'*HC2diff); end
            cnt(3) = cnt(3)+1;
        else
            Hdiff = Hdiff + HC2diff;
        end
    end
    cache.gradFBE = (Hdiff - cache.diff/gam);
end

%% precompute quantities needed for the line search procedures
function [cache, cnt] = CacheLSData(prob, dir, cache)
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];
    cache.dir = dir;
    if prob.istheref1
        if prob.isthereC1
            if prob.isC1fun, cache.C1dir = prob.C1(dir);
            else cache.C1dir = prob.C1*dir; end
            if prob.isQfun, cache.QC1dir = prob.Q(cache.C1dir);
            else cache.QC1dir = prob.Q*cache.C1dir; end
            if prob.isC1fun, cache.C1tQC1dir = prob.C1t(cache.QC1dir);
            else cache.C1tQC1dir = prob.C1'*cache.QC1dir; end
            cnt(2) = cnt(2)+2;
        else
            cache.C1dir = dir;
            if prob.isQfun, cache.QC1dir = prob.Q(cache.C1dir);
            else cache.QC1dir = prob.Q*cache.C1dir; end
            cache.C1tQC1dir = cache.QC1dir;
        end
        cnt(1) = cnt(1)+1;
        cache.f1linear = cache.gradf1x'*dir;
        cache.f1quad = cache.C1dir'*cache.QC1dir;
    end
    if prob.istheref2
        if prob.isthereC2
            if prob.isC2fun, cache.C2dir = prob.C2(dir);
            else cache.C2dir = prob.C2*dir; end
            cnt(3) = cnt(3)+1;
        else
            cache.C2dir = dir;
        end
    end
    if prob.istherelin
        cache.lindir = prob.l'*dir;
    end
end
