classdef FBS < forbes.solvers.AbstractIterativeSolver
    properties
        f1, A1, f2, A2, g, x0
        Lf, gam
        x, z, z_prev
        A1x, gradf1_A1x, f1_A1x
        A2x, gradf2_A2x, f2_A2x
        A1z_prev, gradf1_A1z_prev, A2z_prev
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
            fprintf('%8s | %11s | %11s\n', 'iter', 'gam', 'fpr');
        end
        function display_progress(obj)
            fprintf('%8d | %8.5e | %8.5e\n', obj.it, obj.gam, norm(obj.x - obj.z, inf));
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
