% QUADRATIC Quadratic function

classdef Quadratic < forbes.functions.Proximable
    properties
        Q, q % Hessian and linear term
        L_prox % stores Cholesky factor for prox
        gam_prox % stores most recently used parameter for prox
        L_conj
        p_conj
        flag_sparse
        flag_large
    end
    methods
        function obj = Quadratic(Q, q)
            if issparse(Q)
                obj.flag_sparse = true;
            else
                obj.flag_sparse = false;
            end
            if size(Q, 2) > 5000
                obj.flag_large = true;
            else
                obj.flag_large = false;
            end
            obj.Q = Q;
            obj.q = q;
            obj.gam_prox = 0;
        end
        function p = is_quadratic(obj)
            p = true;
        end
        function p = is_convex(obj)
            % TODO implement this for real
            p = true;
        end
        function p = is_strongly_convex(obj)
            % TODO implement this for real
            p = true;
        end
    end
end
