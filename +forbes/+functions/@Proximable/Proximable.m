classdef Proximable < handle
    properties
        cnt_prox
        cnt_gradient
    end
    methods
        function obj = Proximable()
            obj.cnt_prox = 0;
            obj.cnt_gradient = 0;
        end
        function reset(obj)
            obj.cnt_prox = 0;
            obj.cnt_gradient = 0;
        end
        function v = call(obj, x)
            [~, v] = obj.gradient(x);
        end
        function [p, v] = prox(obj, x, gamma)
            [p, v] = obj.compute_prox(x, gamma);
            obj.cnt_prox = obj.cnt_prox + 1;
        end
        function [g, v] = gradient(obj, x)
            [g, v] = obj.compute_gradient(x);
            obj.cnt_gradient = obj.cnt_gradient + 1;
        end
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
        function p = is_null(obj)
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
