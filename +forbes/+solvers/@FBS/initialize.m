function initialize(obj, f, A, g, x0)
    obj.f = f;
    obj.A = A;
    obj.g = g;
    obj.x0 = x0;

    % compute (approximate) Lipschitz constant (if necessary)

    if ~isinf(obj.opt.Lf)
        if isfield(obj.opt, 'adaptive')
            obj.adaptive = obj.opt.adaptive;
        else
            obj.adaptive = false;
        end
        % nothing to compute if Lipschitz constant is provided by the user
        obj.Lf = obj.opt.Lf;
    else
        if isfield(obj.opt, 'adaptive')
            obj.adaptive = obj.opt.adaptive;
        else
            obj.adaptive = ~f.is_quadratic();
        end
        if ~f.is_quadratic() || obj.adaptive
            obj.Lf = forbes.utils.lipschitz_uppbnd(f, A, x0);
        else
            obj.Lf = forbes.utils.lipschitz_quadratic(f, A, x0);
        end
    end

    % set stepsize, initialize vectors

    obj.gam = 1.0/obj.Lf;
    obj.x = obj.x0;
    obj.v = obj.x0;
end
