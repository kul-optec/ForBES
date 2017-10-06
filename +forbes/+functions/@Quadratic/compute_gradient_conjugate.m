function [g, v] = compute_gradient_conjugate(obj, y)
    if isempty(obj.L_conj)
        factor_gradient_conjugate(obj);
    end
    if obj.flag_sparse
        rhs = y-obj.q;
        g(obj.p_conj,1) = obj.L_conj'\(obj.L_conj\rhs(obj.p_conj));
        v = 0.5*(y-obj.q)'*g;
    else
        g = obj.L_conj'\(obj.L_conj\(y-obj.q));
        v = 0.5*(y-obj.q)'*g;
    end
end

function factor_gradient_conjugate(obj)
    if obj.flag_sparse
        [obj.L_conj,flag,obj.p_conj] = chol(obj.Q,'lower','vector');
        if flag~=0
            error('Q is not positive definite')
        end
    else
        [obj.L_conj,flag] = chol(obj.Q,'lower');
        if flag~=0
            error('Q is not positive definite')
        end
    end
end
