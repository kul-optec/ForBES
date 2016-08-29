%QUADRATICOVERAFFINE Allocates a quadratic function over an affine subspace.
%
%   QUADRATICOVERAFFINE(A, b, Q, q) builds the function
%
%       f(x) = 0.5*(x'*Q*x)+q'*x subject to A*x = b

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

function obj = quadraticOverAffine(A, b, Q, q)
    if nargin < 2, error('you should provide at least 2 arguments'); end
    if nargin < 3, Q = 1; end
    if nargin < 4, q = zeros(size(A, 2)); end
    if isscalar(Q), Q = Q*speye(size(A, 2)); end
    obj.isQuadratic = 0;
    obj.isConjQuadratic = 1;
    obj.makefconj = @() make_quadraticOverAffine_conj(Q, q, A, b);
end


function fc = make_quadraticOverAffine_conj(Q, q, A, b)
    m = size(A,1);
    [L,D,P] = ldl([Q A';A sparse(m,m)]);
    fc = @(y) eval_quadraticOverAffine_conj(y, L, D, P, Q, q, b);
end

function [val, grad] = eval_quadraticOverAffine_conj(y, L, D, P, Q, q, b)
    grad = P*(L'\(D\(L\(P'*[y-q; b]))));
    grad = grad(1:length(q),1);
    val = -(0.5*grad'*(Q*grad)+(q-y)'*grad);
end
