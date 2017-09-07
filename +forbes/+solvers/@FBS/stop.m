function flag_stop = stop(obj)
    if norm(obj.FPR_y, inf)/obj.gam <= obj.opt.tol
        flag_stop = true;
    else
        flag_stop = false;
    end
end
