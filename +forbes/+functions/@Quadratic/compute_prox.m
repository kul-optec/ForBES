function [p, v] = compute_prox(obj, x, gam)
    n = length(x);
    I = speye(n);
    if obj.flag_large
        [p, flag] = pcg(I + gam*obj.Q, x - gam*obj.q);
    else
        if gam ~= obj.gam_prox % more robust test?
            % factor matrix when gam changes (or at first call)
            obj.gam_prox = gam;
            if obj.flag_sparse
                obj.L_prox = ldlchol(I + gam*obj.Q);
            else
                obj.L_prox = chol(I + gam*obj.Q); % do differently for sparse?
            end
        end
        if obj.flag_sparse
            p = ldlsolve(obj.L_prox, x - gam*obj.q);
        else
            p = obj.L_prox\(obj.L_prox'\(x - gam*obj.q));
        end
    end
    v = 0.5*(p'*(obj.Q*p)) + obj.q'*p; % can we save something here?
end
