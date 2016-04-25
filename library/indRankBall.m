%INDRANKBALL Allocates the nuclear norm function
%
%   INDRANKBALL(m, n, r) builds the function
%       
%       g(X) = 0 if rank(X) <= r, +infinity otherwise
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

function obj = indRankBall(m, n, r, method)
    if nargin < 3
        error('you must provide the number of rows and columns, m and n, and rank r as arguments');
    end
    if nargin < 4, method = 'svds'; end
    switch method
        case 'svds'
            obj.makeprox = @() @(x, gam) call_indRankBall_proj_svds(x, m, n, r);
        case 'lansvd'
            obj.makeprox = @() @(x, gam) call_indRankBall_proj_lansvd(x, m, n, r);
        otherwise
            error('unknown method for computing SVDs');
    end
    obj.isConvex = 0;
end

function [prox, val] = call_indRankBall_proj_svds(x, m, n, r)
    [U, S, V] = svds(reshape(x, m, n), r, 'L');
    prox = reshape(U*(S*V'), m*n, 1);
    if nargout >= 2
        val = 0;
    end
end

function [prox, val] = call_indRankBall_proj_lansvd(x, m, n, r)
    [U, S, V] = lansvd(reshape(x, m, n), r, 'L');
    prox = reshape(U*(S*V'), m*n, 1);
    if nargout >= 2
        val = 0;
    end
end
