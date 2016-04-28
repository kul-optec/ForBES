%QUADLOSSOVERAFFINE Allocates the squared distance function over an affine subspace.
%
%   QUADLOSSOVERAFFINE(A, b, w, p) returns the function
%       
%       f(x) = 0.5*sum_i w_i(x_i-p_i)^2 subject to A*x = b
%
%   Requires LDLCHOL and LDLSOLVE from SuiteSparse by Tim Davis.
%   See: http://faculty.cse.tamu.edu/davis/suitesparse.html
%
% Copyright (C) 2015, Lorenzo Stella and Panagiotis Patrinos
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

function obj = quadLossOverAffine(A, b, w, p)
    if nargin < 4, p = 0; end
    if nargin < 3, w = 1; end
    if isscalar(w), w = repmat(w, size(A,2), 1); end
    obj.isConjQuadratic = 1;
    obj.makefconj = @() make_quadLossOverAffine_conj(w, p, A, b);
end

function fun = make_quadLossOverAffine_conj(w, p, A, b)
    LD = ldlchol(A*diag(sparse(1./sqrt(w))), 1e-12);
    fun = @(y) call_quadLossOverAffine_conj(y, LD, w, p, A, b);
end

function [val, grad] = call_quadLossOverAffine_conj(y, LD, w, p, A, b)
    wyp = (y./w)+p;
    grad = wyp-(A'*ldlsolve(LD, A*wyp-b))./w;
    gradp = grad-p;
    val = y'*grad-0.5*(gradp'*(w.*gradp));
end

