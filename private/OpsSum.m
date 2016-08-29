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

function ops = OpsSum(ops1, ops2)

ops.C1      = ops1.C1 + ops2.C1;
ops.C2      = ops1.C2 + ops2.C2;
ops.f1      = ops1.f1 + ops2.f1;
ops.gradf1  = ops1.gradf1 + ops2.gradf1;
ops.f2      = ops1.f2 + ops2.f2;
ops.gradf2  = ops1.gradf2 + ops2.gradf2;
ops.g       = ops1.g + ops2.g;
ops.proxg   = ops1.proxg + ops2.proxg;
