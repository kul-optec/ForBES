function [p, v] = prox(obj, x, gam)
    lamgam = obj.lam*gam;
    if gam ~= obj.gam_prox % more robust test?
        % factor matrix when gam changes (or at first call)
        obj.gam_prox = gam;
        I = speye(size(obj.S));
        t0 = tic();
        if obj.flag_sparse
            if size(obj.A, 1) <= size(obj.A, 2)
                obj.L_prox = ldlchol(obj.A, 1/lamgam);
            else
                obj.L_prox = ldlchol(obj.A', 1/lamgam);
            end
        else
            obj.L_prox = chol(obj.S + I/lamgam);
        end
        obj.time_fact = toc(t0);
    end
    q = obj.Atb + x/lamgam;
    if size(obj.A, 1) <= size(obj.A, 2)
        if obj.flag_sparse
            s = ldlsolve(obj.L_prox, obj.A*q);
        else
            s = obj.L_prox\(obj.L_prox'\(obj.A*q));
        end
        p = lamgam*(q - obj.A'*s);
    else
        if obj.flag_sparse
            p = ldlsolve(obj.L_prox, q);
        else
            p = obj.L_prox\(obj.L_prox'\q);
        end
    end
    v = (obj.lam/2)*norm(obj.A*p-obj.b, 2)^2;
end
