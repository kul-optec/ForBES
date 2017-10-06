% INDBOX Indicator function of a box

classdef IndBox < forbes.functions.Proximable
    properties
        lo, hi % box boundaries
    end
    methods
        function obj = IndBox(lo, hi)
            obj.lo = lo;
            obj.hi = hi;
        end
    end
end
