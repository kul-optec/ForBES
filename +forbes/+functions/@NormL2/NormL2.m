% NORML2 L2 (Euclidean) norm

classdef NormL2 < forbes.functions.Proximable
    properties
        lam % coefficient
    end
    methods
        function obj = NormL2(lam)
            if nargin < 1, lam = 1.0; end
            obj.lam = lam;
        end
        function p = is_convex(obj)
            p = true;
        end
    end
end
