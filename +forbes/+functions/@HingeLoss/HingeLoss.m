% HINGELOSS Hinge loss function

classdef HingeLoss < forbes.functions.Proximable
    properties
        mu, b
    end
    methods
        function obj = HingeLoss(mu, b)
            if nargin < 1, mu = 1.0; end
            if nargin < 2, b = 1.0; end
            obj.mu = mu;
            obj.b = b;
        end
        function p = is_convex(obj)
            p = true;
        end
    end
end
