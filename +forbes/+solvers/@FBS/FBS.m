classdef FBS < forbes.solvers.AbstractIterativeSolver
    properties
        f, A, g, x0
        Lf, gam
        x, v, y, z
        FPR_y
        opt, adaptive
    end
    methods
        function obj = FBS(varargin)
            opt = struct(varargin{:});
            opt = forbes.solvers.FBS.defaults(opt);
            obj@forbes.solvers.AbstractIterativeSolver(opt.maxit, opt.verbose);
            obj.opt = opt;
        end
        function display_header(obj)
            fprintf('%8s | %8s\n', 'iter', 'fpr');
        end
        function display_progress(obj)
            fprintf('%8d | %8.5e\n', obj.it, norm(obj.FPR_y, inf));
        end
    end
    methods (Static)
        function opt = defaults(opt)
            default_opt.verbose = false;
            default_opt.maxit = 10000;
            default_opt.tol = 1e-5;
            default_opt.Lf = inf;
            default_opt.fast = false;
            default_fields = fieldnames(default_opt);
            for i = 1:length(default_fields)
                k = default_fields{i};
                v = getfield(default_opt, k);
                if ~isfield(opt, k), opt = setfield(opt, k, v); end
            end
        end
    end
end
