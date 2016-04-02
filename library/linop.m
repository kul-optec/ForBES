%LINOP Allocates a linear operator given function handles computing the
%mapping and adjoint mapping.
%
%   LINOP(op, adj, m, n) builds the linear operator defined by MATLAB
%   function handles op (computing the operator) and adj (computing the
%   adjoint operator). Parameters m and n are respectively the dimensions
%   of the target space and domain of op.
%   
%   All parameters are compulsory.
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

function obj = linOp(op, adj, m, n)
    if nargin < 4
        error('you should provide 4 arguments: op, adj, m, n');
    end
    if ~isa(op, 'function_handle') || ~isa(adj, 'function_handle')
        error('first two arguments should be function handles');
    end
    obj.m = m;
    obj.n = n;
    obj.makeop = @() op;
    obj.makeadj = @() adj;
end
