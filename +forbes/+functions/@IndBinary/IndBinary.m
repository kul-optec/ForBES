% INDBINARY Indicator function the set {low,high}^n 

classdef IndBinary < forbes.functions.Proximable
    properties
        low, high
    end
    methods
        function obj = IndBinary(low, high)
            if nargin == 0
                low = 0;
                high = 1;
            end
            obj.low = low;
            obj.high = high;
        end
    end
end
