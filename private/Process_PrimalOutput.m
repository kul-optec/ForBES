function out1 = Process_PrimalOutput(prob, out)

out1.name = out.name;
out1.message = out.message;
out1.flag = out.flag;
out1.gam = out.gam;
[out1.x1, out1.x2, out1.z] = prob.Get_DualPoints(out.x, out.gam);
out1.y = out.x;
out1.iterations = out.iterations;
if isfield(out, 'operations'), out1.operations = out.operations; end
if isfield(out, 'record'), out1.record = out.record; end
out1.residual = out.residual;
out1.ts = out.ts;
out1.dual = out;
out1.prob = prob;

end
