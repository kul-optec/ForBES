% INDRANKBALL Indicator function of the set of matrices with rank at most r

classdef IndRankBall < forbes.functions.Proximable
    properties
        r % ball radius
        method % indicates which svd to use
    end
    methods
        function obj = IndRankBall(r, method)
            if nargin < 2
                method = 'svds';
            end
            switch method
                case 'svds'
                    obj.method = 1;
                case 'lansvd'
                    obj.method = 2;
                otherwise
                    error('unknown method');
            end
            obj.r = r;
        end
    end
end
