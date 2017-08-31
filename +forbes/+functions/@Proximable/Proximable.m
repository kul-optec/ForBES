classdef Proximable < handle
    methods
        function p = is_prox_accurate(obj)
            p = true;
        end
        function p = is_separable(obj)
            p = false;
        end
        function p = is_convex(obj)
            p = false;
        end
        function p = is_singleton(obj)
            p = false;
        end
        function p = is_cone(obj)
            p = false;
        end
        function p = is_affine(obj)
            p = obj.is_singleton();
        end
        function p = is_set(obj)
            p = obj.is_affine() || obj.is_cone();
        end
        function p = is_smooth(obj)
            p = obj.is_quadratic();
        end
        function p = is_quadratic(obj)
            p = false;
        end
        function p = is_generalized_quadratic(obj)
            p = obj.is_quadratic() || obj.is_affine();
        end
        function p = is_strongly_convex(obj)
            p = false;
        end
        function p = has_hessian(obj)
            p = false;
        end
    end
    methods (Static)
        function mu = get_gram_diagonal(M)
            mus = M*(M'*ones(size(M, 1), 1));
            if (max(mus)-min(mus))/(1+abs(min(mus))) > 1e-14
                mu = 0;
            else
                mu = 0.5*(max(mus)+min(mus));
            end
        end
    end
end
