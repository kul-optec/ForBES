% FORBES_QP
%
%   FORBES_QP(H, q, A, lb, ub, Aeq, beq, lx, ux, opt, out1) solves the
%   quadratic programming problem
%
%     minimize   (1/2)*x'*H*x + q'*x
%     subject to     Aeq*x == beq
%                lb <= A*x <= ub
%                lx <=   x <= ux
%
%   All arguments are optional and may be empty, except for the first.
%   If last argument out1 is specified, then the solution process is
%   warm-started based on the output of a previous call.

function out = forbes_qp(H, q, A, lb, ub, Aeq, beq, lx, ux, opt, out1)

    t0 = tic();

    % Arguments parsing

    if nargin < 1 || isempty(H)
        error('first parameter H is mandatory');
    end
    
    n = size(H, 2);
    
    if nargin < 2 || isempty(q)
        q = zeros(n, 1);
    end
    
    if nargin < 3 || isempty(A)
        flag_ineq = 0;
    else
        flag_ineq = 1;
        if size(A, 2) ~= n
            error('argument A is incompatible with H and q');
        end
    end
    
    m = size(A, 1);
    
    if nargin < 4 || isempty(lb)
        lb = -inf(m, 1);
    end
    
    if any(size(lb) ~= [m, 1])
        error('argument lb is incompatible with A');
    end
    
    if nargin < 5 || isempty(ub)
        ub = +inf(m, 1);
    end
    
    if any(size(ub) ~= [m, 1])
        error('argument ub is incompatible with A');
    end
    
    if nargin < 6 || isempty(Aeq)
        flag_eq = 0;
    else
        flag_eq = 1;
        if size(Aeq, 2) ~= n
            error('size of Aeq is incompatible with H and q');
        end
    end
    
    if flag_eq == 1 && (nargin < 7 || isempty(beq) || any(size(beq) ~= [size(Aeq,1), 1]))
        error('argument beq is incompatible with Aeq');
    end

    if nargin < 8 || isempty(lx)
        lx = -inf(n, 1);
    else
        if isscalar(lx)
            lx = lx*ones(n, 1);
        end
    end
    
    if nargin < 9 || isempty(ux)
        ux = +inf(n, 1);
    else
        if isscalar(ux)
            ux = ux*ones(n, 1);
        end
    end
    
    if nargin < 10, opt = []; end
    if nargin < 11, out1 = []; end
    
    if ~isfield(opt, 'prescale') || isempty(opt.prescale)
        opt.prescale = true;
    end
    
    % Problem setup and solution
    
    if flag_ineq == 0 && flag_eq == 0
        f = quadratic(H, q);
        g = indBox(lx, ux);
        if isempty(out1)
            x0 = zeros(n, 1);
        else
            opt.Lf = out1.solver.prob.Lf;
            x0 = out1.x;
        end
        tprep = toc(t0);
        out = forbes(f, g, x0, [], [], opt);
    else
        A_ext = A;
        lb_ext = lb;
        ub_ext = ub;
        % Extend inequality constraints so as to include bounds on x
        % (This should only be done if necessary)
        A_ext = [A_ext; speye(n)];
        lb_ext = [lb_ext; lx];
        ub_ext = [ub_ext; ux];
        m_ext = m + n;
        % Scale inequality constraints
        if opt.prescale
%             scaling_A = 1./sum(A.^2, 2);
%             scaling_A = 1./max(abs(A),[],2);
            scale = 1./sqrt(diag(A_ext*(H\A_ext')));
            A_ext = diag(sparse(scale))*A_ext;
            lb_ext = scale.*lb_ext;
            ub_ext = scale.*ub_ext;
        end
        if flag_eq == 0
            f = quadratic(H, q);
            opt_eigs.issym = 1;
            opt_eigs.tol = 1e-3;
            if (n <= 500 && min(eig(H)) <= 1e-16) || ...
               (n >  500 && eigs(H, 1, 'SM', opt_eigs) <= 1e-16)
                out.status = 2;
                out.msg = 'not strongly convex';
                return;
            end
        else
            f = quadraticOverAffine(Aeq, beq, H, q);
        end
        g = indBox(lb_ext, ub_ext);
        constr = {A_ext, -speye(m_ext), zeros(m_ext, 1)};
        if isempty(out1)
            % cold start
            y0 = zeros(m_ext, 1);
        else
            % warm start
            opt.Lf = out1.solver.dual.prob.Lf;
            y0 = [out1.y_ineq; out1.y_bnd];
        end
        tprep = toc(t0);
        out_forbes = forbes(f, g, y0, [], constr, opt);
    end
    
    ttot = toc(t0);
    
    out.status = out_forbes.flag;
    out.msg = out_forbes.message;
    out.x = out_forbes.x1;
    out.y_ineq = out_forbes.y(1:m);
    out.y_bnd = out_forbes.y(m+1:end);
    out.pobj = (out.x'*(H*out.x))/2 + q'*out.x;
    out.dobj = -out_forbes.dual.objective(end); % dual is solved as minimization
    out.iterations = out_forbes.iterations;
    out.preprocess = tprep;
    out.time = ttot;
    out.solver = out_forbes;
end