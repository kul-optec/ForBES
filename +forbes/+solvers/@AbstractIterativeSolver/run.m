function run(obj, varargin)
    t0 = tic();
    obj.initialize(varargin{:});
    obj.it = 0;
    obj.status = 1;
    obj.message = 'exceeded max iter';
    obj.inittime = toc(t0);
    if obj.verbose
        obj.display_header();
    end
    while obj.it < obj.maxit
        stop = obj.iterate();
        obj.it = obj.it+1;
        if obj.verbose
            obj.display_progress();
        end
        if stop
            obj.status = 0;
            obj.message = 'converged';
            break;
        end
    end
    obj.runtime = toc(t0);
end
