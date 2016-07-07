function obj = ProcessFunction(obj)
    if ~isfield(obj, 'isConvex') || isempty(obj.isConvex), obj.isConvex = 0; end
    if ~isfield(obj, 'isQuadratic') || isempty(obj.isQuadratic), obj.isQuadratic = 0; end
    if ~isfield(obj, 'isConjQuadratic') || isempty(obj.isConjQuadratic), obj.isConjQuadratic = 0; end
    if ~isfield(obj, 'hasHessian') || isempty(obj.hasHessian), obj.hasHessian = 0; end
end