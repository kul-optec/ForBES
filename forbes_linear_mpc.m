% FORBES_LINEAR_MPC
%
%   FORBES_LINEAR_MPC(A, B, Q, R, Q_N, g_x, L_x, g_u, L_u, opt)

function out = forbes_linear_mpc(x0, Q, R, Q_N, A, B, N, g_x, L_x, g_u, L_u, opt, out_prev)

    t0 = tic();

    % Arguments parsing

    if nargin < 7
      error('some mandatory arguments are missing');
    end

    n_x = size(Q, 1);
    n_u = size(R, 1);

    if ~exist('g_x', 'var') || isempty(g_x), g_x = zeroFunction(); end
    if ~exist('L_x', 'var') || isempty(L_x), L_x = speye(n_x); end
    if ~exist('g_u', 'var') || isempty(g_u), g_u = zeroFunction(); end
    if ~exist('L_u', 'var') || isempty(L_u), L_u = speye(n_u); end

    m_x = size(L_x, 1);
    m_u = size(L_u, 1);
    
    % Compute total dimensions
    
    n_x_tot = (N+1)*n_x;
    n_u_tot = N*n_u;
    m_x_tot = (N+1)*m_x;
    m_u_tot = N*m_u;

    % Problem setup and solution (multiple shooting)

    f = lqrCost(x0, Q, R, Q_N, A, B, N);
    g = separableSum({g_x, g_u}, {m_x_tot, m_u_tot});
    
    % Build big constraint matrix
    
    L = sparse(m_x_tot+m_u_tot, n_x_tot+n_u_tot);
    for i = 0:N-1
        base_i = i*m_x; base_j = i*(n_x+n_u);
        L(base_i+(1:m_x), base_j+(1:n_x)) = L_x;
        base_i = m_x_tot+i*m_u; base_j = i*(n_x+n_u)+n_x;
        L(base_i+(1:m_u), base_j+(1:n_u)) = L_u;
    end
    base_i = N*m_x; base_j = N*(n_x+n_u);
    L(base_i+(1:m_x), base_j+(1:n_x)) = L_x;
    
    % Now the problem to solve is
    % 
    %   minimize f(xu) + g(z) subject to L*xu = z
    
    if ~exist('out_prev', 'var') || isempty(out_prev)
        y0 = zeros(m_x_tot+m_u_tot, 1);
    else
        y0 = out_prev.y;
    end
    
    tpre = toc(t0);
    
    out_forbes = forbes(f, g, y0, [], {L, -1, zeros(m_x_tot+m_u_tot,1)}, opt);

    ttot = toc(t0);
    
    out.xu = out_forbes.x1;
    out.z = out_forbes.z;
    out.y = out_forbes.y;
    out.solver = out_forbes;
    out.preprocess = tpre;
    out.time = ttot;

end

