% DISTL2 L2 (Euclidean) distance from a convex set

classdef DistL2 < forbes.functions.Proximable
    properties
        ind % indicator function of the (convex) set
        lam
    end
    methods
        function obj = DistL2(ind, lam)
            if nargin < 2, lam = 1.0; end
            obj.ind = ind;
            obj.lam = lam;
        end
        function p = is_convex(obj)
            p = true;
        end
    end
end
