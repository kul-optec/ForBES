%MATFAC
%
%   MATFAC(A, r) returns the function
%
%       f(x) = 0.5*||A-UV||^2_F
%
%   where x = [U(:); V(:)]. If A is n-times-m then U is n-times-r and V is
%   r-times-m, therefore it must be length(x) = (n+m)*r.

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

function obj = matFac(A, r)
    if nargin < 2
        error('two arguments are required: A (matrix to factor) and r (rank of the factorization)');
    end
    [n, m] = size(A);
    obj.makef = @() @(x) call_matrixFactorization_fun(x, A, n, r, m);
end

function [val, grad] = call_matrixFactorization_fun(x, A, n, r, m)
    nr = n*r;
    mr = m*r;
    U = reshape(x(1:nr), n, r);
    V = reshape(x(nr+1:nr+mr), r, m);
    res = U*V - A;
    val = 0.5*norm(res, 'fro')^2;
    grad = [reshape((res*V'), nr, 1); reshape((U'*res), mr, 1)];
end
