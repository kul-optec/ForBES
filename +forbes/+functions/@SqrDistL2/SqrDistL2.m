% SQRDISTL2 Squared L2 (Euclidean) distance from a convex set

classdef SqrDistL2 < forbes.functions.Proximable
    properties
        ind % indicator function of the (convex) set
    end
    methods
        function obj = SqrDistL2(ind)
            obj.ind = ind;
        end
    end
end
