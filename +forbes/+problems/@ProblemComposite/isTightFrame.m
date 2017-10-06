function [b, mu] = isTightFrame(prob, L)
  if isa(L, 'function')
    % how to handle this case?
    b = false;
    mu = 0;
  elseif isa(L, 'numeric')
    diagLLt = sum(L.*L, 2);
    if (max(diagLLt)-min(diagLLt))/max(diagLLt) > 10*eps || diagLLt(1) <= 10*eps
      b = false;
      mu = 0;
    else
      b = true;
      mu = diagLLt(1);
    end
end
