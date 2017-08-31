% INDBALLL2 Indicator function of the L2 (Euclidean) ball

classdef IndBallL2 < forbes.functions.Proximable
    properties
        R % ball radius
    end
    methods
        function obj = IndBallL2(R)
            if nargin < 1
                R = 1.0;
            end
            obj.R = R;
        end
    end
end
