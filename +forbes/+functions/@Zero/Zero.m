% ZERO Null function

classdef Zero < forbes.functions.Proximable
    methods
        function obj = Zero()
        end
        function p = is_convex(obj)
            p = true;
        end
        function p = is_quadratic(obj)
            p = true;
        end
    end
end
