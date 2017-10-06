% HUBERLOSS Huber loss function

classdef HuberLoss < forbes.functions.Proximable
    properties
        del % coefficient
    end
    methods
        function obj = HuberLoss(del)
            obj.del = del;
        end
        function p = is_convex(obj)
            p = true;
        end
        function p = is_smooth(obj)
            p = true;
        end
    end
end
