% SCALETRANSLATE Scales and translates a function input
%
%   SCALETRANSLATE(f, a, b) returns the function f(a.*x - b).
%
%   Only scalar a (i.e., uniform scaling) is currently available.

classdef ScaleTranslate < forbes.functions.Proximable
    properties
        f, a, b
    end
    methods
        function obj = ScaleTranslate(f, a, b)
            if nargin < 1, a = 1.0; end
            if nargin < 2, b = 0.0; end
            if ~isscalar(a)
                error('only uniform scaling is currently available; a should be a scalar');
            end
            % TODO: check that a ~= 0.0
            obj.f = f;
            obj.a = a;
            obj.b = b;
        end
        function p = is_generalized_quadratic(obj)
            p = obj.f.is_generalized_quadratic();
        end
        function p = is_quadratic(obj)
            p = obj.f.is_quadratic();
        end
        function p = is_smooth(obj)
            p = obj.f.is_smooth();
        end
        function p = is_convex(obj)
            p = obj.f.is_convex();
        end
        function p = is_strongly_convex(obj)
            p = obj.f.is_strongly_convex();
        end
        function p = is_singleton(obj)
            p = obj.f.is_singleton();
        end
        function p = is_cone(obj)
            p = obj.f.is_cone();
        end
        function p = is_affine(obj)
            p = obj.f.is_affine();
        end
        function p = is_set(obj)
            p = obj.f.is_set();
        end
    end
end
