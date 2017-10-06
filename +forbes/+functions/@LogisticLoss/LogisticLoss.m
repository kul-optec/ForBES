% LOGISTICLOSS Logistic loss function

classdef LogisticLoss < forbes.functions.Proximable
    properties
        mu % coefficient
    end
    methods
        function obj = LogisticLoss(mu)
            if nargin < 1, mu = 1; end
            obj.mu = mu;
        end
        function p = is_convex(obj)
            p = true;
        end
        function p = is_smooth(obj)
            p = true;
        end
        function p = has_hessian(obj)
            p = true;
        end
    end
end
