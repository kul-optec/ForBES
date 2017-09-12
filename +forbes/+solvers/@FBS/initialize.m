function initialize(obj, f1, A1, f2, A2, g, x0)
    obj.f1 = f1;
    obj.A1 = A1;
    obj.f2 = f2;
    obj.A2 = A2;
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
            obj.adaptive = ~f2.is_null();
        end
        if ~f2.is_null() || obj.adaptive
            f = forbes.functions.SeparableSum({f1, f2}, {size(A1, 1), size(A2, 1)});
            A = [A1; A2];
            obj.Lf = forbes.utils.lipschitz_lowbnd(f, A, x0);
        else
            obj.Lf = forbes.utils.lipschitz_quadratic(f1, A1, x0);
        end
    end

    % set stepsize, initialize vectors

    obj.gam = 1.0/obj.Lf;
    obj.x = obj.x0;
end
