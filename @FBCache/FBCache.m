classdef FBCache < handle
  properties
    % cache parameters
    prob, x, gam
    % intermediate arrays
    res1x, f1x, gradf1res1x, gradf1x
    res2x, f2x, gradf2res2x, gradf2x
    gradfx
    y, z, FPR
    dir1, dir2
    C1dir1, QC1dir1, C1tQC1dir1, C2dir1
    C1dir2, QC1dir2, C1tQC1dir2, C2dir2
    gradFBE
    % scalars
    fx, gz, normFPR, FBE, dFBE
    flinx
    f1linear1, f1quad1, lindir1
    f1linear2, f1quad2, lindir2
    f1cross
    % flags
    flagEvalf, flagGradStep, flagProxGradStep, flagFBE, flagGradFBE
    flagLineSearch1, flagLineSearch2
    % operation counter
    ops, flagOps
  end
  methods
    function cache = FBCache(myprob, myx, mygam, myops)
      cache.prob = myprob;
      cache.x = myx;
      cache.gam = mygam;
      cache.flagEvalf = false;
      cache.flagGradStep = false;
      cache.flagProxGradStep = false;
      cache.flagFBE = false;
      cache.flagGradFBE = false;
      cache.flagLineSearch1 = false;
      cache.flagLineSearch2 = false;
      if nargin < 4 || isempty(myops)
          cache.ops = [];
          cache.flagOps = false;
      else
          cache.ops = myops;
          cache.flagOps = true;
      end
    end
  end
end
