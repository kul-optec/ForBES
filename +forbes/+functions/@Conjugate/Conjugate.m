% CONJUGATE Fenchel conjugate function (of some other given function)

classdef Conjugate < forbes.functions.Proximable
    properties
        f
    end
    methods
        function obj = Conjugate(f)
            obj.f = f;
        end
        function p = is_convex(obj)
            p = true;
        end
        function p = is_strongly_convex(obj)
            p = obj.f.is_smooth();
        end
        function p = is_smooth(obj)
            p = obj.f.is_strongly_convex();
        end
        function p = is_quadratic(obj)
            p = obj.f.is_strongly_convex() && obj.f.is_generalized_quadratic();
        end
        function p = is_generalized_quadratic(obj)
            p = obj.f.is_quadratic();
        end
    end
end
