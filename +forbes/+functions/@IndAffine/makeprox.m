function makeprox(obj)
    disp('invoked makeprox!');
    % can we make this more efficient, e.g., for sparse matrices?
    obj.L_prox = chol(obj.A*obj.A','lower');
end
