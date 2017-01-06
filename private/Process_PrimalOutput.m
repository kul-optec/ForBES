function out = Process_PrimalOutput(prob, dualprob, dualout)
  
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
out.dual = dualout;
out.dual.prob = dualprob;

end
