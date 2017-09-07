function u = mtimes(obj, v)
  if obj.currmem > 0
    H = obj.YS(obj.idx)/(obj.Y(:,obj.idx)'*obj.Y(:,obj.idx));
    u = lbfgs_mex(obj.S, obj.Y, obj.YS, H, v, int32(obj.idx), int32(obj.currmem));
  else
    u = v;
  end
end
