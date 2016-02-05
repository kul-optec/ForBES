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

% initialize output stuff
ts = zeros(1, opt.maxit);
residual = zeros(1, opt.maxit);
msgTerm = '';

% initialize operations counter
ops = OpsInit();

% initialize stuff
gam = SelectGamma(prob, opt);

% display header
if opt.display >= 2
    fprintf('%6s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'tau');
end

alpha = 0.5;
c = 0.9;
thetabar = 0.1;
flag = -1; % to track what happened at every iteration
flag_gamma = 0;

MAXIMUM_Lf = 1e15;
% MINIMUM_tau = 1e-15;
MINIMUM_d = 1e-15;

t0 = tic();

cache_current = CacheInit(prob, prob.x0, gam);
[cache_current, ops1] = CacheProxGradStep(cache_current, gam);
ops = OpsSum(ops, ops1);
cache_previous = [];

eta = cache_current.normdiff;

for it = 1:opt.maxit
    
    flag_gamma = 0;
    
    % backtracking on gamma
    
    cache_z = CacheInit(prob, cache_current.z, gam);
    [cache_z, ops1] = CacheEvalf(cache_z);
    ops = OpsSum(ops, ops1);
    
    if prob.unknownLf || opt.adaptive
        while cache_z.fx > cache_current.fx + cache_current.gradfx'*cache_current.diff ...
                + prob.Lf/2*cache_current.normdiff^2
            if prob.Lf >= MAXIMUM_Lf, break; end
            prob.Lf = 2*prob.Lf;
            gam = SelectGamma(prob, opt);
            flag_gamma = 1;
            [cache_current, ops1] = CacheProxGradStep(cache_current, gam);
            ops = OpsSum(ops, ops1);
            cache_z = CacheInit(prob, cache_current.z, gam);
            [cache_z, ops1] = CacheEvalf(cache_z);
            ops = OpsSum(ops, ops1);
        end
    end
    
    if prob.Lf >= MAXIMUM_Lf
        msgTerm = ['estimate for Lf became too large: ', num2str(prob.Lf)];
        flagTerm = 1;
        break;
    end
    
    % adjust sigma
    
    sig = SelectSigma(prob, opt, gam);
    
    % trace stuff
    
    ts(1, it) = toc(t0);
    residual(1, it) = norm(cache_current.diff, 'inf');
    
    % check for termination
    
    if isnan(cache_current.normdiff)
        msgTerm = 'something went wrong';
        flagTerm = 1;
        break;
    end
    if ~flag_gamma
        if residual(1, it) <= opt.tol
            msgTerm = 'reached optimum (up to tolerance)';
            flagTerm = 0;
            break;
        end
    end
    
    % select a direction
    
    switch opt.method
        case 11
            if it == 1 || flag_gamma
                B = eye(prob.n);
                d = cache_current.diff;
            else
                s = cache_current.x - cache_previous.x;
                y = cache_previous.diff-cache_current.diff ;
                rho = (B\y)'*s/(s'*s);
                if abs(rho) >= thetabar, theta = 1.0;
                else theta = (1-sign(rho)*thetabar)/(1+rho); end
                B = B + (theta/(s'*s))*((y-B*s)*s');
                d = B\cache_current.diff;
            end
        case 12 % BFGS
            opt.optsL.UT = true; opt.optsL.TRANSA = true;
            opt.optsU.UT = true;
            if it == 1 || flag_gamma
                d = cache_current.diff;
                R = eye(prob.n);
            else
                %%% x' - x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Sk = cache_current.x - cache_previous.x;
                Yk = cache_previous.diff-cache_current.diff ;
                %%% other options (is this additional gradient eval needed?)
                YSk = Yk'*Sk;
                Bs = R'*(R*Sk);
                sBs = Sk'*Bs;
                if YSk>= 0.2*sBs
                    theta = 1;
                else
                    theta = (0.8*sBs)/(sBs-YSk);
                end
                r = theta*Yk + (1-theta)*Bs;
                R = cholupdate(cholupdate(R,r/sqrt(Sk'*r)),Bs/sqrt(sBs),'-');
%                 if YSk > 0
%                     R = cholupdate(cholupdate(R,Yk/sqrt(YSk)),Bs/sqrt(sBs),'-');
%                 else
%                     skipCount = skipCount+1;
%                 end
                d = linsolve(R,linsolve(R,cache_current.diff,opt.optsL),opt.optsU);
                %                     dir = -R\(R'\cache_current.gradFBE);
            end
        case 13
            if it == 1 || flag_gamma
                d = cache_current.diff;
                skipCount = 0;
                LBFGS_col = 1;
                LBFGS_mem = 0;
            else
                %%% x' - x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Sk = cache_current.x - cache_previous.x;
                Yk = cache_previous.diff-cache_current.diff ;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                YSk = Yk'*Sk;
                if YSk > 1e-12;
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
                    d = LBFGS(S, Y, YS, H, cache_current.diff, int32(LBFGS_col), int32(LBFGS_mem));
                else
                    d = cache_current.diff;
                end
            end
        otherwise
            error('search direction not implemented');
    end
    
    % select a stepsize
    
    if norm(d) <= MINIMUM_d
        % if d is zero then tau = 0.0 and the algorithm becomes FBS
        tau = 0.0;
    else
        tau = 1.0;
    end
    
    while 1
        x = cache_current.z + tau*(d - cache_current.diff);
        cache_next = CacheInit(prob, x, gam);
        [cache_next, ops1] = CacheProxGradStep(cache_next, gam);
        ops = OpsSum(ops, ops1);
        if cache_current.normdiff <= c*eta && cache_next.normdiff <= c*eta
            eta = cache_current.normdiff;
            cache_previous = cache_current;
            cache_current = cache_next;
            flag = 1;
            break;
        end
        [cache_current, ops1] = CacheFBE(cache_current, gam);
        ops = OpsSum(ops, ops1);
        [cache_next, ops1] = CacheFBE(cache_next, gam);
        ops = OpsSum(ops, ops1);
        if cache_next.FBE <= cache_current.FBE - sig*cache_current.normdiff^2
            cache_previous = cache_current;
            cache_current = cache_next;
            flag = 0;
            break;
        end
        tau = alpha*tau;
    end

    % display stuff
    if opt.display == 1
        if mod(it, 100) == 0
            fprintf('.');
        end
        if mod(it, 4000) == 0
            fprintf('\n');
        end
    elseif opt.display >= 2
        fprintf('%6d %7.4e %7.4e %7.4e %d\n', it, gam, residual(1,it), tau, flag);
    end

end

if it == opt.maxit
    msgTerm = 'exceeded maximum iterations';
    flagTerm = 1;
end

% pack up results
out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
out.x = cache_current.z;
out.iterations = it;
out.operations = ops;
out.residual = residual(1, 1:it);
out.ts = ts(1, 1:it);
out.prob = prob;
out.opt = opt;
out.gam = gam;

function gam = SelectGamma(prob, opt)

gam = 0.95/prob.Lf;

function sig = SelectSigma(prob, opt, gam)

sig = (1-gam*prob.Lf)/(4*gam);

