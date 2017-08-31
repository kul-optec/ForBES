% SQRNORML2 Squared L2 (Euclidean) norm

classdef SqrNormL2 < forbes.functions.Proximable
    properties
        w % weight(s)
    end
    methods
        function obj = SqrNormL2(w)
            if nargin < 1, w = 1; end
            obj.w = w;
        end
        function p = is_convex(obj)
            p = true;
        end
        function p = is_quadratic(obj)
            p = true;
        end
        function p = has_hessian(obj)
            p = true;
        end
    end
end
