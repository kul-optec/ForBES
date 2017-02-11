% FORBES_LINEAR_MPC
%
%   FORBES_LINEAR_MPC(mpc_prob, opt) solves
%   the linear model predictive control problem
%
%   min. 0.5*sum((x[k]-xref)'*Q*(x[k]-xref) + u[k]'*R*u[k], k=0,...,N-1) [stage cost]
%
%        + 0.5*((x[N]-xref)'*Q_N*(x[N]-xref))                            [final cost]
%
%        + sum(g_s(L_s(x[k], u[k])), k=0,...,N-1)          [stage penalty]
%
%        + g_N(L_N(x[N]))                                  [final penalty]
%
%        + g_c(L_c(x, u)                                   [coupling term]
%
%   s.t. x[0] = x0
%
%        x[k+1] = A x[k] + B u[k], k = 0,...,N-1           [dynamics]
% 
%   mpc_prob is a structure containing the following problem parameters:
% 
%       mpc_prob.x0, mpc_prob.xref, mpc_prob.Q, mpc_prob.R, mpc_prob.Q_N
%       mpc_prob.A, mpc_prob.B, mpc_prob.N, mpc_prob.L_s, mpc_prob.L_N
% 
%   penalty functions g_s, g_N are determined as follows:
% 
%       mpc_prob.s_min, mpc_prob.s_max: lower/upper bound on the stage 
%       mpc_prob.x_N_min, mpc_prob.x_N_max: lower/upper bound on final state 
% 
%   and
% 
%       mpc_prob.stage_w = [w_1, ..., w_{m_s}], the weights to apply
%       to the linear penalty for each constraint violation (+inf: hard
%       constraint)
% 
%       mpc_prob.final_w = [w_1, ..., w_{m_N}], analogous to the
%       previous case
% 

% function out = forbes_linear_mpc(x0, xref, Q, R, Q_N, A, B, N, g, L, opt, out_prev)
function out = forbes_linear_mpc(mpc_prob, opt, out_prev)

    t0 = tic();

    % Arguments parsing

    if nargin < 1
      error('you must provide an mpc_prob structure');
    end

    if ~exist('opt','var'), opt = []; end
    
    if ~isfield(opt,'prescale') || isempty(opt.prescale)
        opt.prescale = 1;
    end

    % Compute problem size

    n_x = size(mpc_prob.A, 2);
    n_u = size(mpc_prob.B, 2);

    m_stage = size(mpc_prob.L_s, 1);
    
    % Make objective term f

    if ~isfield(mpc_prob, 'xref') || isempty(mpc_prob.xref)
        f = lqrCost(mpc_prob.x0, mpc_prob.Q, mpc_prob.R, mpc_prob.Q_N, ...
            mpc_prob.A, mpc_prob.B, mpc_prob.N);
        mpc_prob.xref = zeros(n_x, 1);
    else
        f = lqrCost(mpc_prob.x0, mpc_prob.Q, mpc_prob.R, mpc_prob.Q_N, ...
            mpc_prob.A, mpc_prob.B, mpc_prob.N, mpc_prob.xref);
    end

    % Build big constraint matrix
    
    if isfield(mpc_prob, 'x_N_ellipse')
        mpc_prob.L_N = mpc_prob.x_N_ellipse{1};
        alpha = mpc_prob.x_N_ellipse{2};
    end
    
    m_final = size(mpc_prob.L_N, 1);

    diag_L = {};
    for k = 1:mpc_prob.N
        diag_L{k} = mpc_prob.L_s;
    end
    diag_L{mpc_prob.N+1} = mpc_prob.L_N;
    L = sparse(blkdiag(diag_L{:})); % ugly
    
    % Compute scaling

    if opt.prescale
        scale = zeros(size(L, 1), 1);
        callfconj = f.makefconj();
        [~, p] = callfconj(zeros(size(L, 2),1));
        for i = 1:size(L, 1)
            [~, dgradi] = callfconj(L(i, :)');
            w = L(i, :)*(dgradi-p);
            if w >= 1e-14
                scale(i) = 1/sqrt(w);
            else
                scale(i) = 1;
            end
        end
    else
        scale = ones(size(L, 1), 1);
    end

    % Scale nonsmooth term & equality constraints matrix

    xu_min = []; xu_max = []; w = [];
    for k=0:mpc_prob.N-1
        xu_min = [xu_min; mpc_prob.s_min];
        xu_max = [xu_max; mpc_prob.s_max];
        w = [w; mpc_prob.stage_w];
    end
    
    if isfield(mpc_prob, 'x_N_ellipse')
        % Case where ellipsoidal final constraint is selected
        xu_min_scaled = scale(1:mpc_prob.N*m_stage).*xu_min;
        xu_max_scaled = scale(1:mpc_prob.N*m_stage).*xu_max;
        w_scaled = w./scale(1:mpc_prob.N*m_stage);
        g_s = distBox(xu_min_scaled, xu_max_scaled, w_scaled);
        % do not scale final constraint
        scale_N = mean(scale(end-m_final+1:end));
        scale(end-m_final+1:end) = scale_N;
        g_N = indBall_l2(scale_N*sqrt(2*alpha), scale_N*mpc_prob.xref);
        g = separableSum({g_s, g_N}, {mpc_prob.N*m_stage, n_x});
    else
        % Case where ordinary (soft/hard) final constraint is selected
        xu_min = [xu_min; mpc_prob.x_N_min];
        xu_max = [xu_max; mpc_prob.x_N_max];
        w = [w; mpc_prob.final_w];
        xu_min_scaled = scale.*xu_min;
        xu_max_scaled = scale.*xu_max;
        w_scaled = w./scale;
        g = distBox(xu_min_scaled, xu_max_scaled, w_scaled);
    end
    
    L_scaled = sparse(diag(scale))*L;

    % Now the problem to solve is
    %
    %   minimize f(xu) + g(z) subject to L_scaled * xu = z

    % Set starting (dual) point

    if ~exist('out_prev', 'var') || isempty(out_prev)
        y0 = zeros(size(L_scaled, 1), 1);
    else
        y0 = out_prev.y;
    end

    tpre = toc(t0);

    out_forbes = forbes(f, g, y0, [], {L_scaled, -1, zeros(length(y0), 1)}, opt);

    ttot = toc(t0);

    out.xu = out_forbes.x1;
    temp = reshape(out_forbes.x1(1:end-n_x), n_x+n_u, mpc_prob.N);
    out.x = [temp(1:n_x,:), out_forbes.x1(end-n_x+1:end)];
    out.u = temp(n_x+1:end,:);
    out.z = out_forbes.z;
    out.y = out_forbes.y; % dual variables
    out.forbes = out_forbes;
    out.preprocess = tpre;
    out.time = ttot;
    out.scaling = scale;

end
