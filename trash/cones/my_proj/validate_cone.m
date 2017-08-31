function K = validate_cone(K)

if (~isfield(K, 'zero')), K.zero = 0; end
if (~isfield(K, 'nn')), K.nn = 0; end
if (~isfield(K, 'soc')), K.soc = []; end
if (~isfield(K, 'sdc')), K.sdc = []; end
if (~isfield(K, 'exp')), K.exp = 0; end

end
