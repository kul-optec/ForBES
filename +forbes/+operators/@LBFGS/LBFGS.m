classdef LBFGS < handle
% LBFGS creates a limited-memory BFGS operator
%
%   H = LBFGS(mem) constructs a limited-memory BFGS operator with memory mem.
%
%   H.push(s, y) pushes a new pair of vectors (s, y) into the buffer.
%
%   u = H*v applies H to vector v; use -(H*v) if the minus is needed.
  properties
    mem % memory
    S, Y, YS, H0, idx, currmem
    cntskip
  end
  methods
    function obj = LBFGS(mem)
      obj.mem = mem;
      obj.S = [];
      obj.Y = [];
      obj.YS = [];
      obj.idx = 0;
      obj.currmem = 0;
      obj.cntskip = 0;
    end
    function reset(obj)
      obj.idx = 0;
      obj.currmem = 0;
    end
  end
end
