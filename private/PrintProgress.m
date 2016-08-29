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

function PrintProgress(varargin)
    if length(varargin) == 1
        it = varargin{1};
        if mod(it, 100) == 0
            fprintf('.');
        end
        if mod(it, 4000) == 0
            fprintf('\n');
        end
    end
    if length(varargin) == 2
        flag = varargin{2};
        if flag == 0
            fprintf('*\n');
        else
            fprintf('!!!\n');
        end
    end
end
