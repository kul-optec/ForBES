function [Lf, adaptive] = Get_Lipschitz(prob, opt)

if isfield(opt, 'Lf') && ~isempty(opt.Lf)
  Lf = opt.Lf;
  if isfield(opt, 'adaptive') && ~isempty(opt.adaptive)
    adaptive = opt.adaptive;
  else
    adaptive = false;
  end
else
  [Lf, exactLf] = prob.EstimateLipschitz();
  if isfield(opt, 'adaptive') && ~isempty(opt.adaptive)
    adaptive = opt.adaptive;
  else
    adaptive = ~exactLf;
  end
end
