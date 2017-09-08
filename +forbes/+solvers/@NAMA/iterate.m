function stop = iterate(obj)
    % TODO: optimize
    if obj.it == 0 || ~obj.adaptive
        obj.A1xk = obj.A1*obj.xk;
        [obj.gradf1_A1xk, obj.f1_A1xk] = obj.f1.gradient(obj.A1xk);
        obj.A2xk = obj.A2*obj.xk;
        [obj.gradf2_A2xk, obj.f2_A2xk] = obj.f2.gradient(obj.A2xk);
    end
    At_gradf_Axk = obj.A1'*obj.gradf1_A1xk + obj.A2'*obj.gradf2_A2xk;
    f_Axk = obj.f1_A1xk + obj.f2_A2xk;
    [obj.xbark, g_xbark] = obj.g.prox(obj.xk - obj.gam*At_gradf_Axk, obj.gam);

    obj.FPR_xk = obj.xk - obj.xbark;
    normFPR_xk = norm(obj.FPR_xk, 'fro');

    uppbnd = f_Axk - At_gradf_Axk(:)'*obj.FPR_xk(:) + 0.5/obj.gam*normFPR_xk^2;

    reset = false;

    if obj.adaptive
        f_Axbark = +inf;
        % TODO: avoid infinite loops here
        while f_Axbark > uppbnd
            A1xbark = obj.A1*obj.xbark;
            A2xbark = obj.A2*obj.xbark;
            [gradf1_A1xbark, f1_A1xbark] = obj.f1.gradient(A1xbark);
            [gradf2_A2xbark, f2_A2xbark] = obj.f2.gradient(A2xbark);
            f_Axbark = f1_A1xbark + f2_A2xbark;
            if f_Axbark > uppbnd + 1e-6*abs(f_Axk)
                reset = true;
                obj.Lf = 2*obj.Lf;
                obj.gam = (1-obj.opt.bet)/obj.Lf;
                [obj.xbark, g_xbark] = obj.g.prox(obj.xk - obj.gam*At_gradf_Axk, obj.gam);
                obj.FPR_xk = obj.xk - obj.xbark;
                normFPR_xk = norm(obj.FPR_xk, 'fro');
                uppbnd = f_Axk - At_gradf_Axk(:)'*obj.FPR_xk(:) + 0.5/obj.gam*normFPR_xk^2;
            end
        end
    else
        A1xbark = 0.0;
        A2xbark = 0.0;
    end

    if obj.stop()
        stop = true;
        return;
    else
        stop = false;
    end

    FBE_xk = uppbnd + g_xbark;

    % Compute direction

    if reset == true
        obj.Hk.reset();
    end

    dk = -(obj.Hk*obj.FPR_xk);

    % Perform backtracking: this looks messy, but it's really simple
    % The mess comes from optimizing the calls to A1, A2 and f1.gradient,
    % and from having to test obj.gam (in the adaptive case).

    tau = 1.0;

    xkdk = obj.xk + dk;
    A1xkdk = obj.A1xk + obj.A1*dk;
    A2xkdk = obj.A2xk + obj.A2*dk;
    [gradf1_A1xkdk, f1_A1xkdk] = obj.f1.gradient(A1xkdk);
    A1t_gradf1_A1xkdk = obj.A1'*gradf1_A1xkdk;

    lin_coeff = 0.0;
    quad_coeff = 0.0;
    A1t_gradf1_A1xbark = 0.0;

    for lsit = 1:obj.opt.maxbacktrack
        wk = (1-tau)*obj.xbark + tau*xkdk;
        A2wk = (1-tau)*A2xbark + tau*A2xkdk;
        A1t_gradf1_A1wk = (1-tau)*A1t_gradf1_A1xbark + tau*A1t_gradf1_A1xkdk;

        % Explanation of next (and some of previous) line(s).
        %
        % Function f1 is quadratic, therefore the following expansion is exact:
        %
        %   f1(x+y) = f1(x) + f1'(x)*y + 0.5*f1"*||y||^2.
        %
        % In this case, we want to compute f1 at A1wk, and to do so we can
        % expand f1 around A1xkdk. In fact
        %
        %   A1wk = (1-tau)*A1xbark + tau*A1xkdk
        %        = A1xkdk + (1-tau)*(A1xbark - A1xkdk).
        %
        % Therefore at every backtracking iteration (i.e. for every tau)
        %
        %   f1(A1wk) = f1(A1xkdk) + f1'(A1xkdk)*(1-tau)*(A1xbark - A1xkdk)
        %              + 0.5*f1"*(1-tau)^2*||A1xbark - A1xkdk||^2.
        %
        % How do we compute f1"? Well, for any two points x, y, one has
        %
        %   f1"*(x-y) = f1'(x) - f1'(y).

        f1_A1wk = f1_A1xkdk + lin_coeff*(1-tau) + quad_coeff*(1-tau)^2;
        [gradf2_A2wk, f2_A2wk] = obj.f2.gradient(A2wk);
        A2t_gradf2_A2wk = obj.A2'*gradf2_A2wk;
        At_gradf_Awk = A1t_gradf1_A1wk + A2t_gradf2_A2wk;
        f_Awk = f1_A1wk + f2_A2wk;
        [wbark, g_wbark] = obj.g.prox(wk - obj.gam*At_gradf_Awk, obj.gam);
        FPR_wk = wk - wbark;
        normFPR_wk = norm(FPR_wk, 'fro');
        reset = false;
        uppbnd = f_Awk - At_gradf_Awk(:)'*FPR_wk(:) + 0.5/obj.gam*normFPR_wk^2;

        reset = false;
        if obj.adaptive
            f_Awbark = +inf;
            % TODO: avoid infinite loops here
            while f_Awbark > uppbnd
                A1wbark = obj.A1*wbark;
                A2wbark = obj.A2*wbark;
                [gradf1_A1wbark, f1_A1wbark] = obj.f1.gradient(A1wbark);
                [gradf2_A2wbark, f2_A2wbark] = obj.f2.gradient(A2wbark);
                f_Awbark = f1_A1wbark + f2_A2wbark;
                if f_Awbark > uppbnd + 1e-6*abs(f_Awk)
                    reset = true;
                    obj.Lf = 2*obj.Lf;
                    obj.gam = (1-obj.opt.bet)/obj.Lf;
                    [wbark, g_wbark] = obj.g.prox(obj.wk - obj.gam*At_gradf_Awk, obj.gam);
                    FPR_wk = obj.wk - obj.wbark;
                    normFPR_wk = norm(FPR_wk, 'fro');
                    uppbnd = f_Awk - At_gradf_Awk(:)'*FPR_wk(:) + 0.5/obj.gam*normFPR_wk^2;
                end
            end
        end
        if reset == true, break; end

        FBE_wk = uppbnd + g_wbark;
        if FBE_wk <= FBE_xk
            xk_backup = obj.xk;
            obj.xk = wbark;
            if obj.adaptive
                obj.A1xk = A1wbark;
                obj.gradf1_A1xk = gradf1_A1wbark;
                obj.f1_A1xk = f1_A1wbark;
                obj.A2xk = A2wbark;
                obj.gradf2_A2xk = gradf2_A2wbark;
                obj.f2_A2xk = f2_A2wbark;
            end
            break;
        end
        if lsit == 1
            if ~obj.adaptive % otherwise we have already computed this stuff
                A1xbark = obj.A1*obj.xbark;
                A2xbark = obj.A2*obj.xbark;
                [gradf1_A1xbark, ~] = obj.f1.gradient(A1xbark);
            end
            A1t_gradf1_A1xbark = obj.A1'*gradf1_A1xbark;
            temp1 = A1xbark - A1xkdk;
            lin_coeff = gradf1_A1xkdk(:)'*temp1(:);
            temp2 = gradf1_A1xbark - gradf1_A1xkdk;
            quad_coeff = 0.5*(temp1(:)'*temp2(:));
        end
        if lsit == obj.opt.maxbacktrack
            xk_backup = obj.xk;
            obj.xk = obj.xbark;
            if obj.adaptive
                obj.A1xk = A1xbark;
                obj.gradf1_A1xk = gradf1_A1xbark;
                obj.f1_A1xk = f1_A1xbark;
                obj.A2xk = A2xbark;
                obj.gradf2_A2xk = gradf2_A2xbark;
                obj.f2_A2xk = f2_A2xbark;
            end
            obj.num_lsfails = obj.num_lsfails + 1;
        end
        tau = 0.5*tau;
    end

    if reset == true
        obj.Hk.reset();
    else
        obj.Hk.push(wk - xk_backup, FPR_wk - obj.FPR_xk);
    end
end
