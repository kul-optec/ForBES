classdef NAMA < forbes.solvers.AbstractIterativeSolver
    properties
        f1, A1, f2, A2, g, x0
        Lf, gam
        xk, xbark
        A1xk, gradf1_A1xk, f1_A1xk % these are useful in the adaptive case
        A2xk, gradf2_A2xk, f2_A2xk %
        FPR_xk, Hk
        opt, adaptive
        num_lsfails
    end
    methods
        function obj = NAMA(varargin)
            opt = struct(varargin{:});
            opt = forbes.solvers.NAMA.defaults(opt);
            obj@forbes.solvers.AbstractIterativeSolver(opt.maxit, opt.verbose);
            obj.opt = opt;
            obj.Hk = opt.method;
            obj.num_lsfails = 0;
        end
        function display_header(obj)
            fprintf('%8s | %8s\n', 'iter', 'fpr');
        end
        function display_progress(obj)
            fprintf('%8d | %8.5e\n', obj.it, norm(obj.FPR_xk, inf));
        end
    end
    methods (Static)
        function opt = defaults(opt)
            default_opt.verbose = false;
            default_opt.maxit = 10000;
            default_opt.tol = 1e-5;
            default_opt.bet = 0.05;
            default_opt.Lf = inf;
            default_opt.method = forbes.directions.LBFGS(10);
            default_opt.maxbacktrack = 10;
            default_fields = fieldnames(default_opt);
            for i = 1:length(default_fields)
                k = default_fields{i};
                v = getfield(default_opt, k);
                if ~isfield(opt, k), opt = setfield(opt, k, v); end
            end
        end
    end
end
