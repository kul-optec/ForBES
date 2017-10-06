% NUCLEARNORM Nuclear norm function

classdef NuclearNorm < forbes.functions.Proximable
    properties
        lam, mode, method, nsv
    end
    methods
        function obj = NuclearNorm(lam, mode, method)
            if nargin < 2, mode = 'exact'; end
            if nargin < 3, method = 'svds'; end
            obj.lam = lam;
            obj.mode = mode;
            obj.method = method;
            obj.nsv = 10;
        end
    end
end
