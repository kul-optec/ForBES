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

function out = GetPrimalOutput(prob, dualprob, dualout)
    out.name = dualout.name;
    out.message = dualout.message;
    out.flag = dualout.flag;
    out.gam = dualout.gam;
    Ax = 0;
    y = dualout.x;
    if isfield(prob, 'f1')
        if isa(prob.A1, 'function_handle')
            [~, out.x1] = dualprob.callf1(-prob.A1t(y));
            Ax = Ax+prob.A1(out.x1);
        else
            [~, out.x1] = dualprob.callf1(-prob.A1'*y);
            Ax = Ax+prob.A1*out.x1;
        end
    end
    if isfield(prob, 'f2')
        if isa(prob.A2, 'function_handle')
            [~, out.x2] = dualprob.callf2(-prob.A2t(y));
            Ax = Ax+prob.A2(out.x2);
        else
            [~, out.x2] = dualprob.callf2(-prob.A2'*y);
            Ax = Ax+prob.A2*out.x2;
        end
    end
    w = y+out.gam*(Ax-prob.b);
    [out.z, ~] = prob.callg(-(prob.B'*w)/out.gam, 1/(prob.muB*out.gam));
    out.y = y;
    out.iterations = dualout.iterations;
    if isfield(dualout, 'operations'), out.operations = dualout.operations; end
    if isfield(dualout, 'record'), out.record = dualout.record; end
    out.residual = dualout.residual;
    out.ts = dualout.ts;
    out.prob = prob;
    out.dual = dualout;
end
