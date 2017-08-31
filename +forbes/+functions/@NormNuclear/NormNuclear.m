% NORMNUCLEAR Nuclear norm function

classdef NormNuclear < forbes.functions.Proximable
    properties
        lam, mode, method
    end
    methods
        function obj = NormNuclear(lam, mode, method)
            if nargin < 2, mode = 'exact'; end
            if nargin < 3, method = 'svds'; end
            obj.lam = lam;
            obj.mode = mode;
            obj.method = method;
        end
    end
end