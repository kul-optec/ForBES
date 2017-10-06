% PRECOMPOSE Precomposition with a linear operator
%
%   PRECOMPOSE(f, L, mu) returns the function f(L*x), where L is a linear
%   operator such that L*L' = mu*I, for mu > 0.

classdef Precompose < forbes.functions.Proximable
    properties
        f, L, mu
    end
    methods
        function obj = Precompose(f, L, mu)
% Here L is a Gram-diagonal operator
% i.e. such that L*L' = mu*I, mu > 0
            if isscalar(L) && nargin > 2
                warning('when L is scalar, mu = L^2; ignoring third argument');
                mu = L^2;
            end
            obj.f = f;
            obj.L = L;
            obj.mu = mu;
        end
    end
end
