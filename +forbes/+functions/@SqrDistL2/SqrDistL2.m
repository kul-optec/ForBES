% SQRDISTL2 Squared L2 (Euclidean) distance from a convex set

classdef SqrDistL2 < forbes.functions.Proximable
    properties
        ind % indicator function of the (convex) set
        lam
    end
    methods
        function obj = SqrDistL2(ind, lam)
            if nargin < 2, lam = 1.0; end
            obj.ind = ind;
            obj.lam = lam;
        end
        function p = is_convex(obj)
            p = true;
        end
        function p = is_smooth(obj)
            p = true;
        end
    end
end
