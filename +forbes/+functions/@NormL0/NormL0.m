% NORML0 L0 pseudo-norm

classdef NormL0 < forbes.functions.Proximable
    properties
        lam % coefficient
    end
    methods
        function obj = NormL0(lam)
            if nargin < 1, lam = 1.0; end
            obj.lam = lam;
        end
    end
end
