% SEPARABLESUM Separable sum of functions

classdef SeparableSum < forbes.functions.Proximable
    properties
        fs
        dims
        idx
        dimsum
    end
    methods
        function obj = SeparableSum(fs, dims, idx)
            l = length(fs);
            if nargin < 3
                idx = 1:l;
            end
            for i = 1:length(dims)
                if numel(dims{i}) == 1, dims{i} = [dims{i}, 1]; end
            end
            dimsum = zeros(length(idx), 1);
            dimsum(1) = prod(dims{1});
            for i = 2:length(idx)
                dimsum(i) = dimsum(i-1) + prod(dims{i});
            end
            obj.fs = fs;
            obj.dims = dims;
            obj.idx = idx;
            obj.dimsum = dimsum;
        end
        function p = is_quadratic(obj)
            p = true;
            for i = 1:length(obj.fs)
                if ~obj.fs{i}.is_quadratic()
                    p = false;
                    break;
                end
            end
        end
        function p = is_convex(obj)
            p = true;
            for i = 1:length(obj.fs)
                if ~obj.fs{i}.is_convex()
                    p = false;
                    break;
                end
            end
        end
        function p = is_strongly_convex(obj)
            p = true;
            for i = 1:length(obj.fs)
                if ~obj.fs{i}.is_strongly_convex()
                    p = false;
                    break;
                end
            end
        end
        function p = is_generalized_quadratic(obj)
            p = true;
            for i = 1:length(obj.fs)
                if ~obj.fs{i}.is_generalized_quadratic()
                    p = false;
                    break;
                end
            end
        end
        function p = is_smooth(obj)
            p = true;
            for i = 1:length(obj.fs)
                if ~obj.fs{i}.is_smooth()
                    p = false;
                    break;
                end
            end
        end
    end
end
