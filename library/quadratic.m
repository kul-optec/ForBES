%QUADRATIC Allocates a quadratic function.
%
%   QUADRATIC(Q, q) builds the function
%
%       f(x) = 0.5*x'*Q*x+q'*x
%
%   Both arguments are required. Matrix Q can be a scalar, in which case
%   it is intended to be a diagonal matrix with uniform diagonal elements.

% Copyright (C) 2015-2016, Lorenzo Stella and Panagiotis Patrinos
%
% This file is part of ForBES.
%
% ForBES is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ForBES is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with ForBES. If not, see <http://www.gnu.org/licenses/>.

function obj = quadratic(Q, q)
    obj.isQuadratic = 1;
    obj.isConjQuadratic = 1;
    obj.hasHessian = 1;
    if isa(Q, 'function_handle')
        obj.makef = @() @(x) call_quadratic_fun_handle(Q, q, x);
        obj.Q = Q;
    else
        obj.makef = @() @(x) call_quadratic_fun_matrix(Q, q, x);
        obj.Q = @(x) Q*x;
    end
    obj.makefconj = @() make_quadratic_conj(Q, q);
end

function [v, g, Q] = call_quadratic_fun_matrix(Q, q, x)
    g = Q*x+q;
    v = 0.5*(g+q)'*x;
end

function [v, g, Q] = call_quadratic_fun_handle(Q, q, x)
    g = Q(x)+q;
    v = 0.5*(g+q)'*x;
end

function fun = make_quadratic_conj(Q, q)
    if issparse(Q)
        [L,flag,p] = chol(Q,'lower','vector');
        if flag~=0
            error('Q is not positive definite')
        end
        fun = @(y) call_quadratic_sparse_conj(L, p, q, y);
    else
        [L,flag] = chol(Q,'lower');
        if flag~=0
            error('Q is not positive definite')
        end
        fun = @(y) call_quadratic_dense_conj(L, q, y);
    end
end

function [v, g] = call_quadratic_sparse_conj(L, p, q, y)
    rhs = y-q;
    g(p,1) = L'\(L\rhs(p));
    v = 0.5*(y-q)'*g;
end

function [v, g] = call_quadratic_dense_conj(L, q, y)
    g = L'\(L\(y-q));
    v = 0.5*(y-q)'*g;
end

