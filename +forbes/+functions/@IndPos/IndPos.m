% INDPOS Indicator function of the nonnegative orthant

classdef IndPos < forbes.functions.Proximable
    properties
        lo % lower boundary
    end
    methods
        function obj = IndPos(lo)
            if nargin < 1, lo = 0.0; end
            obj.lo = lo;
        end
    end
end
