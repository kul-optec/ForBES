% in minfbe:    v = gradient of FBE
% in zerofpr:   v = residual
function [dir, cacheDir] = ComputeDir(prob, opt, it, hasGammaChanged, sk, yk, v, cacheDir)

% compute direction according to the method
switch opt.methodID
    case {1, 7} % STEEPEST DESCENT and BARZILAI-BORWEIN
        dir = -v;
    case 2 % BFGS
        if it == 1 || hasGammaChanged
            dir = -v;
            cacheDir.R = eye(prob.n);
        else
            R = cacheDir.R;
            YSk = yk'*sk;
            Bs = R'*(R*sk);
            sBs = sk'*Bs;
            if YSk > 0
                R = cholupdate(cholupdate(R,yk/sqrt(YSk)),Bs/sqrt(sBs),'-');
%             else
%                 cacheDir.cntSkip = cacheDir.cntSkip+1;
            end
            dir = -linsolve(R,linsolve(R,v,opt.optsL),opt.optsU);
            cacheDir.R = R;
        end
    case 3 % L-BFGS
        if it == 1 || hasGammaChanged
            dir = -v; % use steepest descent direction initially
            cacheDir.LBFGS_col = 0; % last column of Sk, Yk that was filled in
            cacheDir.LBFGS_mem = 0; % current memory of the method
        else
            YSk = yk'*sk;
            if YSk > 0
                cacheDir.LBFGS_col = 1+mod(cacheDir.LBFGS_col, opt.memory);
                cacheDir.LBFGS_mem = min(cacheDir.LBFGS_mem+1, opt.memory);
                cacheDir.S(:,cacheDir.LBFGS_col) = sk;
                cacheDir.Y(:,cacheDir.LBFGS_col) = yk;
                cacheDir.YS(cacheDir.LBFGS_col) = YSk;
            else
                cacheDir.cntSkip = cacheDir.cntSkip+1;
            end
            if cacheDir.LBFGS_mem > 0
                H = cacheDir.YS(cacheDir.LBFGS_col)/...
                    (cacheDir.Y(:,cacheDir.LBFGS_col)'*cacheDir.Y(:,cacheDir.LBFGS_col));
                dir = LBFGS(cacheDir.S, cacheDir.Y, cacheDir.YS, H, ...
                    -v, int32(cacheDir.LBFGS_col), int32(cacheDir.LBFGS_mem));
            else
                dir = -v;
            end
        end
    case 8 % Broyden
        if it == 1 || hasGammaChanged % || mod(it, 5) == 1
            dir = -v;
            cacheDir.R = eye(prob.n);
        else
            R = cacheDir.R;
            sts = sk'*sk;
            if opt.bopt == 1 % enforces nonzero determinant
                sig = 0.1;
                prev_v = cacheDir.prev_v;
                prev_tau = cacheDir.prev_tau;
                gammak = ((R*yk)'*sk)/sts;
                if abs(gammak) < sig,
                    theta = (1-sgn(gammak)*sig)/(1-gammak);
                    yk = theta*yk - (1-theta)*prev_tau*prev_v;
                end
            elseif opt.bopt == 2 % enforces positive curvature along sk
                sig = 0.5;
                prev_v = cacheDir.prev_v;
                prev_tau = cacheDir.prev_tau;
                sty = sk'*yk;
                stv = sk'*prev_v;
                if sty < sig*prev_tau*abs(stv)
                    theta = (1+sgn(stv)*sig)*prev_tau*stv/(prev_tau*stv + sty);
                    yk = theta*yk - (1-theta)*prev_tau*prev_v;
                end
            end
            Ry = R*yk;
            R = R + (sk-Ry)*(sk'*R)/(sk'*Ry);
            dir = -R*v;
            cacheDir.R = R;
        end
    case 9 % limited-memory Broyden (using compact representation)
        if it == 1 || hasGammaChanged
            dir = -v; % use steepest descent direction initially
            cacheDir.S = []; % stores vectors sk
            cacheDir.Y = []; % stores vectors yk
            cacheDir.W = []; % stores columns of H0*Y-S
            cacheDir.StY = []; % stores inner products <sk,yk>
            cacheDir.M = []; % these two are m-by-m where m is the memory of the method
            cacheDir.LBroyden_mem = 0;
        else
            % damping
            if opt.bopt == 2 % enforces positive curvature along sk
                sig = 0.1;
                prev_v = cacheDir.prev_v;
                prev_tau = cacheDir.prev_tau;
                sty = sk'*yk;
                stv = sk'*prev_v;
                if sty < sig*prev_tau*abs(stv)
                    theta = (1+sign(stv)*sig)*prev_tau*stv/(prev_tau*stv + sty);
                    yk = theta*yk - (1-theta)*prev_tau*prev_v;
                end
            end
            delta = 1; % diagonal of H0
            wk = delta*yk - sk;
            if cacheDir.LBroyden_mem == opt.memory, idx0 = 2;
            else idx0 = 1; cacheDir.LBroyden_mem = cacheDir.LBroyden_mem+1; end
            S0 = cacheDir.S(:,idx0:end);
            Y0 = cacheDir.Y(:,idx0:end);
            W0 = cacheDir.W(:,idx0:end);
            StY0 = cacheDir.StY(idx0:end,idx0:end);
            M0 = cacheDir.M(idx0:end,idx0:end);
            % update matrices S, Y, W, StY, M
            cacheDir.S = [S0, sk];
            cacheDir.Y = [Y0, yk];
            cacheDir.W = [W0, wk];
            if isempty(Y0), cacheDir.StY = sk'*yk;
            else cacheDir.StY = [[StY0; sk'*Y0], cacheDir.S'*yk]; end
            if isempty(S0), cacheDir.M = 0;
            else cacheDir.M = [[M0; S0(:,end)'*S0], zeros(cacheDir.LBroyden_mem, 1)]; end
            K = delta * cacheDir.StY - cacheDir.M;
            % compute direction
            dir = delta*(cacheDir.W * (K\(cacheDir.S'*v)) - v);
        end
%     case 10 % limited-memory Broyden (recursive formula)
%         if it == 1 || hasGammaChanged
%             dir = -v;
%             cacheDir.LB_col = 0;
%             cacheDir.LB_mem = 0;
%             cacheDir.S = zeros(length(dir), opt.memory);
%             cacheDir.Y = cacheDir.S;
%             cacheDir.HY = cacheDir.S;
%             cacheDir.YS = zeros(opt.memory, 1);
%             cacheDir.SHY = zeros(opt.memory, 1);
%         else
%             % damping
%             if opt.bopt == 2 % enforces positive curvature along sk
%                 sig = 0.1;
%                 prev_v = cacheDir.prev_v;
%                 prev_tau = cacheDir.prev_tau;
%                 sty = sk'*yk;
%                 stv = sk'*prev_v;
%                 if sty < sig*prev_tau*abs(stv)
%                     theta = (1+sign(stv)*sig)*prev_tau*stv/(prev_tau*stv + sty);
%                     yk = theta*yk - (1-theta)*prev_tau*prev_v;
%                 end
%             end
%             cacheDir.LB_col      = 1 + mod(cacheDir.LB_col, opt.memory);
%             cacheDir.LB_mem      = min(cacheDir.LB_mem+1, opt.memory);
%             cacheDir.S(:,cacheDir.LB_col) = sk;
%             cacheDir.Y(:,cacheDir.LB_col) = yk;
%             HYk = yk;
%             for jm = 0:cacheDir.LB_mem-2
%                 j    = mod(cacheDir.LB_col+jm, cacheDir.LB_mem) + 1;
%                 HYj  = cacheDir.HY(:,j);
%                 Sj   = cacheDir.S(:,j);
%                 SHyj = cacheDir.SHY(j);
%                 HYk  = HYk + (Sj-HYj)*((Sj'*HYk)/SHyj);
%             end
%             cacheDir.HY(:,cacheDir.LB_col) = HYk;
%             cacheDir.SHY(cacheDir.LB_col)  = sk'*HYk;
%             dir = -v;
%             for jm = 1:cacheDir.LB_mem-1
%                 j    = mod(cacheDir.LB_col+jm, cacheDir.LB_mem) + 1;
%                 HYj  = cacheDir.HY(:,j);
%                 Sj   = cacheDir.S(:,j);
%                 SHyj = cacheDir.SHY(j);
%                 dir  = dir + (Sj-HYj)*((Sj'*dir)/SHyj);
%             end
%         end
    case 10 % limited-memory Broyden (recursive formula)
        if it == 1 || hasGammaChanged || mod(it, opt.memory) == 1
            dir = -v;
            cacheDir.LBroyden_mem = 0;
            cacheDir.S = [];
            cacheDir.Y = [];
            cacheDir.W = [];
        else
            % damping
            if opt.bopt == 2 % enforces positive curvature along sk
                sig = 0.1;
                prev_v = cacheDir.prev_v;
                prev_tau = cacheDir.prev_tau;
                sty = sk'*yk;
                stv = sk'*prev_v;
                if sty < sig*prev_tau*abs(stv)
                    theta = (1+sign(stv)*sig)*prev_tau*stv/(prev_tau*stv + sty);
                    yk = theta*yk - (1-theta)*prev_tau*prev_v;
                end
            end
            delta = 1; % diagonal of H0
            if cacheDir.LBroyden_mem == opt.memory, idx0 = 2;
            else idx0 = 1; cacheDir.LBroyden_mem = cacheDir.LBroyden_mem+1; end
            S0 = cacheDir.S(:,idx0:end);
            Y0 = cacheDir.Y(:,idx0:end);
            W0 = cacheDir.W(:,idx0:end);
            w = delta*yk;
            dir = -delta*v;
            for j = 1:size(W0, 2)
                w = w + (S0(:,j)'*w)*W0(:,j);
                dir = dir + (S0(:,j)'*dir)*W0(:,j);
            end
            wk = (sk-w)/(sk'*w);
            dir = dir + (sk'*dir)*wk;
            % update matrices S, Y, W
            cacheDir.S = [S0, sk];
            cacheDir.Y = [Y0, yk];
            cacheDir.W = [W0, wk];
        end
    otherwise
        error('search direction not implemented');
end

function s = sgn(a)
s = 2*(a>=0)-1;
