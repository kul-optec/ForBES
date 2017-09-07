function stop = iterate(obj)
    if obj.opt.fast
        theta = 2/(obj.it+2); % since it starts from 0
        obj.y = (1-theta)*obj.x + theta*obj.v;
    else
        obj.y = obj.x;
    end

    [gradf_Ay, f_Ay] = obj.f.gradient(obj.A*obj.y);
    At_gradf_Ay = obj.A'*gradf_Ay;
    [obj.z, ~] = obj.g.prox(obj.y - obj.gam*At_gradf_Ay, obj.gam);

    obj.FPR_y = obj.y - obj.z;
    normFPR_y = norm(obj.FPR_y, 'fro');

    uppbnd = f_Ay - At_gradf_Ay(:)'*obj.FPR_y(:) + 0.5/obj.gam*normFPR_y^2;

    if obj.adaptive
        f_Az = +inf;
        % TODO: avoid infinite loops here
        while f_Az > uppbnd
            [~, f_Az] = obj.f.gradient(obj.A*obj.z);
            if f_Az > uppbnd + 1e-6*abs(f_Ay)
                obj.Lf = 2*obj.Lf;
                obj.gam = 1/obj.Lf;
                [obj.z, ~] = obj.g.prox(obj.y - obj.gam*At_gradf_Ay, obj.gam);
                obj.FPR_y = obj.y - obj.z;
                normFPR_y = norm(obj.FPR_y, 'fro');
                uppbnd = f_Ay - At_gradf_Ay(:)'*obj.FPR_y(:) + 0.5/obj.gam*normFPR_y^2;
            end
        end
    end

    if obj.stop()
        stop = true;
        return;
    else
        stop = false;
    end

    if obj.opt.fast
        obj.v = obj.x + (obj.z-obj.x)/theta;
    end

    obj.x = obj.z;
end
