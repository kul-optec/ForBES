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

function tau0 = ComputeTau0(prob, opt, it, hasGammaChanged, sk, yk, dir)

switch opt.methodID
    case {2, 3, 8, 9, 10} % (limited-memory) quasi-Newton methods
        tau0 = 1.0;
    case 7 % Barzilai-Borwein
        if it == 1 || hasGammaChanged
            tau0 = 1.0/norm(dir, inf);
        else
            tau0 = (sk'*sk)/(sk'*yk);
        end
end
