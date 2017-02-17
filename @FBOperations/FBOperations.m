% FBOPERATIONS
%
%   Class holding counters that keep track of how many operations of the
%   various kinds are executed. For example, FBCache uses an FBOPERATIONS
%   object when evaluating proximal gradient operations, the FBE and so on.

classdef FBOperations < handle
  properties
    C1, f1, gradf1
    C2, f2, gradf2
    g, proxg
  end
  methods
    function ops = FBOperations()
      ops.C1 = 0;
      ops.f1 = 0;
      ops.gradf1 = 0;
      ops.C2 = 0;
      ops.f2 = 0;
      ops.gradf2 = 0;
      ops.g = 0;
      ops.proxg = 0;
    end
    function addC1(ops)
      ops.C1 = ops.C1+1;
    end
    function addf1(ops)
      ops.f1 = ops.f1+1;
    end
    function addgradf1(ops)
      ops.gradf1 = ops.gradf1+1;
    end
    function addC2(ops)
      ops.C2 = ops.C2+1;
    end
    function addf2(ops)
      ops.f2 = ops.f2+1;
    end
    function addgradf2(ops)
      ops.gradf2 = ops.gradf2+1;
    end
    function addg(ops)
      ops.g = ops.g+1;
    end
    function addproxg(ops)
      ops.proxg = ops.proxg+1;
    end
  end
end
