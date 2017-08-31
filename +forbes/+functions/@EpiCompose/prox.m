function [p, v] = prox(obj, s, gam)
    switch obj.flag
        case 1 % f = Quadratic
            if gam ~= obj.gam_prox % more robust test?
                % factor matrix when gam changes (or at first call)
                obj.gam_prox = gam;
                obj.L_prox = chol(obj.Q + (1/gam)*obj.A'*obj.A); % do differently for sparse?
            end
            obj.x = obj.L_prox\(obj.L_prox'\((obj.A'*s)/gam - obj.q));
            v = 0.5*(obj.x'*(obj.Q*obj.x)) + obj.q'*obj.x; % can we save something here?
        case 2 % f = QuadraticOverAffine
            error('not implemented');
        case 3 % f = Any proximable function and A'*A = mu*Id, mu > 0
            [obj.x, v] = obj.f.prox(obj.A'*s/obj.mu, gam/obj.mu);
        otherwise
            error('not implemented');
    end
    p = obj.A*obj.x;
end
