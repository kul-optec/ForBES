% EPICOMPOSE Epi-composition
%
%   EPICOMPOSE(f, A) returns the function g(y) = inf {f(x) : Ax = y}.

classdef EpiCompose < forbes.functions.Proximable
    properties
        f, A
        mu
        flag % f is 1: Quadratic, 2: QuadraticOverAffine, 3: Anything else
        Q, q, C, d % to store data of (generalized) quadratic functions
        gam_prox, L_prox % to store Cholesky factor
        x
    end
    methods
        function obj = EpiCompose(f, A)
            if isa(f, 'Quadratic')
                obj.flag = 1;
                obj.Q = f.Q;
                obj.q = f.q;
                obj.gam_prox = 0;
                if isscalar(A), A = A*speye(size(f.Q, 2)); end
            else
                obj.flag = 3;
                obj.mu = Proximable.get_gram_diagonal(A');
                if obj.mu == 0
                    error('A must be s.t. A''*A = mu*Id');
                end
            end
            obj.f = f;
            obj.A = A;
        end
    end
end
