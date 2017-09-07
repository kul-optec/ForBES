function initialize(obj, f1, A1, f2, A2, g, x0)
    obj.f1 = f1;
    obj.A1 = A1;
    obj.f2 = f2;
    obj.A2 = A2;
    obj.g = g;
    obj.x0 = x0;

    % compute (approximate) Lipschitz constant (if necessary)

    if ~isinf(obj.opt.Lf)
        % nothing to compute if Lipschitz constant is provided by the user
        obj.Lf = obj.opt.Lf;
        if isfield(obj.opt, 'adaptive')
            obj.adaptive = obj.opt.adaptive;
        else
            obj.adaptive = false;
        end
    else
        % TODO: finish here
        error('not implemented');
    end

    % set stepsize, initialize vectors

    obj.gam = (1-obj.opt.bet)/obj.Lf;
    obj.xk = obj.x0;
end
