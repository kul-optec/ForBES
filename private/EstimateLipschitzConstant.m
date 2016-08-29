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

function L = EstimateLipschitzConstant(prob)

    delta = max(1e-12, prob.x0*1e-6);

    if prob.istheref1
        if prob.isthereC1
            C1x0 = prob.C1*prob.x0;
            C1x1 = prob.C1*(prob.x0+delta);
            [~, gradf1C1x0] = prob.callf1(C1x0);
            [~, gradf1C1x1] = prob.callf1(C1x1);
            grad1x0 = prob.C1'*gradf1C1x0;
            grad1x1 = prob.C1'*gradf1C1x1;
        else
            [~, grad1x0] = prob.callf1(prob.x0);
            [~, grad1x1] = prob.callf1(prob.x0+delta);
        end
    else
        grad1x0 = 0;
        grad1x1 = 0;
    end

    if prob.istheref2
        if prob.isthereC2
            C2x0 = prob.C2*prob.x0;
            C2x1 = prob.C2*(prob.x0+delta);
            [~, gradf2C2x0] = prob.callf2(C2x0);
            [~, gradf2C2x1] = prob.callf2(C2x1);
            grad2x0 = prob.C2'*gradf2C2x0;
            grad2x1 = prob.C2'*gradf2C2x1;
        else
            [~, grad2x0] = prob.callf2(prob.x0);
            [~, grad2x1] = prob.callf2(prob.x0+delta);
        end
    else
        grad2x0 = 0;
        grad2x1 = 0;
    end

    L = norm(grad1x1 + grad2x1 - grad1x0 - grad2x0)/norm(delta);

end
