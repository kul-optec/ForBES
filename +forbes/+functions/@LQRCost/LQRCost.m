% LQRCOST LQR cost function

classdef LQRCost < forbes.functions.Proximable
    properties
        x0, Q, R, Q_f, A, B, N, QR
        LRs, Ss, Ks, Ms, Ls
        tilt, diff
    end
    methods
        function obj = LQRCost(x0, Q, R, Q_f, A, B, N, xref)
            obj.x0 = x0;
            obj.Q = Q;
            obj.R = R;
            obj.Q_f = Q_f;
            obj.A = A;
            obj.B = B;
            obj.N = N;
            obj.QR = blkdiag(obj.Q, obj.R);
            if nargin > 7
                obj.tilt = [repmat([obj.Q*xref; zeros(size(obj.R, 1), 1)], obj.N, 1); obj.Q_f*xref];
                obj.diff = (obj.N+1)/2*norm(xref)^2;
            else
                obj.tilt = 0;
                obj.diff = 0;
            end
            obj.RiccatiFactor();
            if exist('RiccatiSolve') ~= 3
                RiccatiPath = fileparts(mfilename('fullpath'));
                mex('-outdir', RiccatiPath, [RiccatiPath, '/RiccatiSolve.c'])
            end
        end
        function RiccatiFactor(obj)
            n = size(obj.Q,1);
            m = size(obj.R,1);
            Ps = zeros(n, n, obj.N+1);
            Ps(:,:,obj.N+1) = obj.Q_f;
            obj.LRs = zeros(m, m, obj.N);
            obj.Ss = zeros(m, n, obj.N);
            obj.Ks = zeros(m, n, obj.N);
            obj.Ms = zeros(m, n, obj.N);
            obj.Ls = zeros(n, n, obj.N);
            for k = obj.N:-1:1
                Rbar = obj.R + obj.B'*(Ps(:,:,k+1)*obj.B);
                Rbar = (Rbar+Rbar')/2;
                LR = chol(Rbar, 'lower');
                obj.LRs(:,:,k) = LR;
                obj.Ss(:,:,k) = obj.B'*(Ps(:,:,k+1)*obj.A);
                obj.Ks(:,:,k) = -(LR'\(LR\obj.Ss(:,:,k)));
                Ps(:,:,k) = obj.Q + obj.A'*(Ps(:,:,k+1)*obj.A) + obj.Ss(:,:,k)'*obj.Ks(:,:,k);
                Ps(:,:,k) = (Ps(:,:,k) + Ps(:,:,k)')/2;
            end
            for k = 1:obj.N
                LR = obj.LRs(:,:,k);
                obj.Ms(:,:,k) = -(LR'\(LR\obj.B'));
                obj.Ls(:,:,k) = (obj.A + obj.B*obj.Ks(:,:,k))';
            end
        end
        function set_x0(obj, x0)
            obj.x0 = x0;
        end
        function p = is_strongly_convex(obj)
            p = true;
        end
        function p = is_generalized_quadratic(obj)
            p = true;
        end
    end
end
