% INDPOINT Indicator function of a singleton

classdef IndPoint < forbes.functions.Proximable
    properties
        p
    end
    methods
        function obj = IndPoint(p)
            obj.p = p;
        end
    end
end
