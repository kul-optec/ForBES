% LEASTSQUARES Construct a linear least-squares cost function

classdef LeastSquares < forbes.functions.Proximable
    properties
        A, b, lam
        Atb
        S % to store square matrix, i.e., A'A or AA' depending on size
        L_prox % to store Cholesky factor for prox computation
        gam_prox
        flag_sparse
        time_fact
    end
    methods
        function obj = LeastSquares(A, b, lam)
            if nargin < 2, b = 0; end
            if nargin < 3, lam = 1.0; end
            obj.A = A;
            obj.b = b;
            obj.lam = lam;
            obj.Atb = A'*b;
            if issparse(A)
                obj.flag_sparse = true;
            else
                obj.flag_sparse = false;
                if size(A, 1) <= size(A, 2)
                    obj.S = A*A';
                else
                    obj.S = A'*A;
                end
            end
            obj.gam_prox = 0;
            obj.time_fact = 0;
        end
        function p = is_convex(obj)
            p = true;
        end
        function p = is_quadratic(obj)
            p = true;
        end
    end
end
