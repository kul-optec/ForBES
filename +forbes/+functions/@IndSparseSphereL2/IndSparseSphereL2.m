% INDSPARSESPHEREL2 Indicator function of the L2 (Euclidean) sphere
% intersected with the L0 pseudo-norm ball of a given radius

classdef IndSparseSphereL2 < forbes.functions.Proximable
    properties
        N % L0 ball radius
        R % L2 sphere radius
    end
    methods
        function obj = IndSparseSphereL2(N, R)
            if nargin < 2
                R = 1.0;
            end
            obj.N = N;
            obj.R = R;
        end
    end
end
