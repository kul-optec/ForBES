%VARIATIONMAP2D Allocates a linear operator computing the discrete gradient
%of a 2D signal (e.g. for image processing applications)
%
%   VARIATIONMAP2D(h, w) returns the linear operator computing the discrete
%   gradient for signals with height h and width w. The signal must be fed
%   to the operator in the stacked-columns form. The output of the operator
%   has size 2hw, and each block of size 2 inside it contains the vertical
%   and horizontal difference relative to the same signal element.
%
%   For elements in the bottom row (rightmost column) the vertical (horizontal)
%   difference is assumed to be zero.
%
%   Example: X is an m-by-n matrix containing grayscale values of an image,
%   then
%   
%       D = variationMap2D(m, n);
%       Dx = D*X(:); % this stacks the columns of X
%
%   and vector Dx (of size 2mn) will contain the finite differences of the
%   elements of X.
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

function obj = variationMap2D(h, w)
    if nargin < 2
        error('you should provide 2 arguments: w (width of the image) and h (height of the image)');
    end
    obj = make_variationMap2D(h, w);
end

function D = make_variationMap2D(h, w)
    D = sparse(2*h*w, h*w);
    blockv = [ones(h-1,1); 0];
    blockh = ones(h,1);
    diagv = kron(ones(w,1), blockv);
    diagh = kron([ones(w-1,1); 0], blockh);
    Dv = spdiags([-diagv, diagv], [0, -1], h*w, h*w)';
    Dh = spdiags([-diagh, diagh], [0, -h], h*w, h*w)';
    D(1:2:end,:) = Dv;
    D(2:2:end,:) = Dh;
end
