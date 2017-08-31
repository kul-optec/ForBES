% NORMLHALF L(1/2) pseudo-norm

classdef NormLhalf < forbes.functions.Proximable
    properties
        lam % coefficient
    end
    methods
        function obj = NormLhalf(lam)
            if nargin < 1, lam = 1.0; end
            obj.lam = lam;
        end
    end
end
