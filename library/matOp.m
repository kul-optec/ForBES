%MATOP Allocates a linear operator given a matrix
%
%   MATOP(A) builds the linear operator associated with A. This is just for
%   completeness, as ForBES actually accepts matrices as linear operators,
%   and treats them accordingly.
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

function obj = matOp(A)
    if nargin < 1
        error('you should provide 1 arguments: A');
    end
    if ~ismatrix(A)
        error('first argument must be a matrix');
    end
    obj.m = size(A, 1);
    obj.n = size(A, 2);
    obj.makeop = @() @(x) A*x;
    obj.makeadj = @() @(y) A'*y;
end
