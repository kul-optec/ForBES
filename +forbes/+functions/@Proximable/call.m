function v = call(obj, x)
    [g, v] = obj.gradient(x);
end
