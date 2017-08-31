% INDSPHEREL2 Indicator function of the L2 (Euclidean) sphere

classdef IndSphereL2 < forbes.functions.Proximable
    properties
        R % ball radius
    end
    methods
        function obj = IndSphereL2(R)
            if nargin < 1
                R = 1.0;
            end
            obj.R = R;
        end
    end
end
