function Util_PrintProgress(varargin)

if length(varargin) == 1
    it = varargin{1};
    if mod(it, 100) == 0
        fprintf('.');
    end
    if mod(it, 4000) == 0
        fprintf('\n');
    end
end

if length(varargin) == 2
    flag = varargin{2};
    if flag == 0
        fprintf('*\n');
    else
        fprintf('!!!\n');
    end
end

end
