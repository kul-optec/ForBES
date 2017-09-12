function stop = iterate(obj)
    if obj.it == 0 || ~obj.adaptive
        obj.A1x = obj.A1*obj.x;
        [obj.gradf1_A1x, obj.f1_A1x] = obj.f1.gradient(obj.A1x);
        obj.A2x = obj.A2*obj.x;
    end
    [obj.gradf2_A2x, obj.f2_A2x] = obj.f2.gradient(obj.A2x);
    At_gradf_Ax = obj.A1'*obj.gradf1_A1x + obj.A2'*obj.gradf2_A2x;
    f_Ax = obj.f1_A1x + obj.f2_A2x;
    obj.z_prev = obj.z;
    [obj.z, ~] = obj.g.prox(obj.x - obj.gam*At_gradf_Ax, obj.gam);

    FPR_x = obj.x - obj.z;
    normFPR_x = norm(FPR_x, 'fro');

    uppbnd = f_Ax - At_gradf_Ax(:)'*FPR_x(:) + 0.5/obj.gam*normFPR_x^2;

    if obj.adaptive
        for it_gam = 1:100
            A1z = obj.A1*obj.z;
            [gradf1_A1z, f1_A1z] = obj.f1.gradient(A1z);
            A2z = obj.A2*obj.z;
            [gradf2_A2z, f2_A2z] = obj.f2.gradient(A2z);
            f_Az = f1_A1z + f2_A2z;
            if f_Az > uppbnd + 1e-6*abs(f_Ax)
                obj.Lf = 2*obj.Lf;
                obj.gam = 1/obj.Lf;
                [z, ~] = obj.g.prox(obj.x - obj.gam*At_gradf_Ax, obj.gam);
                FPR_x = obj.x - obj.z;
                normFPR_x = norm(FPR_x, 'fro');
                uppbnd = f_Ax - At_gradf_Ax(:)'*FPR_x(:) + 0.5/obj.gam*normFPR_x^2;
            else
                break;
            end
        end
    end

    if norm(FPR_x, inf)/obj.gam <= obj.opt.tol
        stop = true;
        return;
    else
        stop = false;
    end

    if obj.it == 0 || obj.opt.fast == false
        obj.x = obj.z;
        if obj.adaptive
            obj.A1x = A1z;
            obj.gradf1_A1x = gradf1_A1z;
            obj.f1_A1x = f1_A1z;
            obj.A2x = A2z;
            obj.gradf2_A2x = gradf2_A2z;
            obj.f2_A2x = f2_A2z;
            obj.A1z_prev = A1z;
            obj.gradf1_A1z_prev = gradf1_A1z;
            obj.A2z_prev = A2z;
        end
    else
        extr = obj.it/(obj.it+3);
        obj.x = obj.z + extr*(obj.z - obj.z_prev);
        if obj.adaptive
            % extrapolate other extrapolable quantities
            diff_A1x = extr*(A1z - obj.A1z_prev);
            obj.A1x = A1z + diff_A1x;
            diff_gradf1_A1x = extr*(gradf1_A1z - obj.gradf1_A1z_prev);
            obj.gradf1_A1x = gradf1_A1z + diff_gradf1_A1x;
            obj.f1_A1x = f1_A1z + gradf1_A1z(:)'*diff_A1x + 0.5*(diff_A1x(:)'*diff_gradf1_A1x(:));
            obj.A2x = A2z + extr*(A2z - obj.A2z_prev);
            % store the z-quantities for future extrapolation
            obj.A1z_prev = A1z;
            obj.gradf1_A1z_prev = gradf1_A1z;
            obj.A2z_prev = A2z;
        end
    end
end
