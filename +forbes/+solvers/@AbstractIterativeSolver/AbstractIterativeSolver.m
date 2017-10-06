classdef AbstractIterativeSolver < handle
    properties
        it, maxit, verbose
        inittime, runtime
        status, message
    end
    methods (Abstract)
        initialize(obj, varargin)
        iterate(obj)
        solution(obj)
    end
    methods
        function obj = AbstractIterativeSolver(varargin)
            obj.maxit = varargin{1};
            obj.verbose = varargin{2};
            obj.status = -1;
            obj.message = 'uninitialized';
        end
        function display_header(obj)
            fprintf('iter\n');
        end
        function display_progress(obj)
            fprintf('%d\n', obj.it);
        end
    end
end
