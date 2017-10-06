function push(obj, s, y)
  ys = y'*s;
  if ys > 0.0
    obj.idx = 1+mod(obj.idx, obj.mem);
    obj.currmem = min(obj.currmem+1, obj.mem);
    obj.S(:,obj.idx) = s;
    obj.Y(:,obj.idx) = y;
    obj.YS(obj.idx) = ys;
  else
    obj.cntskip = obj.cntskip + 1;
  end
end
