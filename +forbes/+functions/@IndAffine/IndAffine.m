% INDAFFINE Indicator function of an affine subspace

classdef IndAffine < forbes.functions.Proximable
    properties
        A, b % define the affine subspace
        L_prox % stores the Cholesky factor needed to compute prox
    end
    methods
        function obj = IndAffine(A, b)
            obj.A = A;
            obj.b = b;
            obj.L_prox = [];
        end
    end
end
