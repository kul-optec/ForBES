% DISTBOXL1 L1 distance from a box

classdef DistBoxL1 < forbes.functions.Proximable
    properties
        lb, ub % box bounds
        w
    end
    methods
        function obj = DistBoxL1(lb, ub, w)
            if nargin < 3, w = 1.0; end
            obj.lb = lb;
            obj.ub = ub;
            obj.w = w;
        end
        function p = is_convex(obj)
            p = true;
        end
    end
end
