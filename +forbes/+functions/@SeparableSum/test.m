function test(obj)
    gam = 0.5 + rand();
    x = [];
    y = [];
    v = 0;
    for i = 1:length(obj.idx)
        x1 = randn(obj.dims{i}, 1);
        x = [x; x1];
        [y1, v1] = obj.fs{obj.idx(i)}.prox(x1, gam);
        y = [y; y1];
        v = v + v1;
    end
    [y_test, v_test] = obj.prox(x, gam);
    assert(norm(y_test - y)/(1+norm(y)) <= 1e-12);
    assert(abs(v_test - v)/(1+abs(v)) <= 1e-12);
end