% Function 0.5||A-UV||^2_F
%   where if A is n-times-m
%   then U is n-times-r and V is r-times-m
%
% The vector of variables x is supposed to be partitioned as
%
%   x = (u, v)
%
% where length(u) = n*r, i.e., u is stacked columns of U
% and length(v) = r*m, i.e., v is stacked columns of V

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
    res = A - U*V;
    val = 0.5*norm(res, 'fro')^2;
    grad = [reshape((U'*res)', nr, 1); reshape((res*V')', mr, 1)];
end
