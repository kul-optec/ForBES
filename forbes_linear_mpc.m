% FORBES_LINEAR_MPC
%
%   FORBES_LINEAR_MPC(x0, A, B, Q, R, Q_N, g_x, L_x, g_u, L_u, opt) solves
%   the linear model predictive control problem
%
%   min. 0.5*sum(x[k]'*Q*x[k] + u[k]'*R*u[k], k=0,...,N-1) [stage cost]
%
%        + 0.5*(x[N]'*Q_N*x[N])                            [final cost]
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

function out = forbes_linear_mpc(x0, Q, R, Q_N, A, B, N, g_s, L_s, g_N, L_N, g_c, L_c, opt, out_prev)

    t0 = tic();

    % Arguments parsing

    if nargin < 11
      error('some mandatory arguments are missing');
    end

    if ~exist('opt','var'), opt = []; end

%     if ~isfield(opt, 'mode'), opt.mode = 'ms'; end

    % Compute problem size

    n_x = size(A, 2);
    n_u = size(B, 2);

    m_stage = size(L_s, 1);
    m_final = size(L_N, 1);

    f = lqrCost(x0, Q, R, Q_N, A, B, N);
    blocks_g = {}; dims_g = {};
    for k = 1:N, blocks_g{k} = g_s; dims_g{k} = m_stage; end
    blocks_g{N+1} = g_N; dims_g{N+1} = m_final;

    % Build big constraint matrix

    diag_L = {};
    for k = 1:N, diag_L{k} = L_s; end
    diag_L{N+1} = L_N;
    L = blkdiag(diag_L{:});

    % Add coupling term if present

    if exist('g_c', 'var') && ~isempty(g_c)
        m_coupling = size(L_c, 1);
        blocks_g{N+2} = g_c;
        dims_g{N+2} = m_coupling;
        L = [L; L_c];
    end

    g = separableSum(blocks_g, dims_g);

    % Now the problem to solve is
    %
    %   minimize f(xu) + h(z) subject to big_L * xu = z

    % Set starting (dual) point

    if ~exist('out_prev', 'var') || isempty(out_prev)
        y0 = zeros(size(L, 1), 1);
    else
        y0 = out_prev.y;
    end

    tpre = toc(t0);

    out_forbes = forbes(f, g, y0, [], {L, -1, zeros(length(y0), 1)}, opt);

    ttot = toc(t0);

    out.xu = out_forbes.x1;
    temp = reshape(out_forbes.x1(1:end-n_x), n_x+n_u, N);
    out.x = [temp(1:n_x,:), out_forbes.x1(end-n_x+1:end)];
    out.u = temp(n_x+1:end,:);
    out.z = out_forbes.z; % slack variables
    out.y = out_forbes.y; % dual variables
    out.solver = out_forbes;
    out.preprocess = tpre;
    out.time = ttot;

end
