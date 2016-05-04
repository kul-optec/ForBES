function outs = forbes_compare(fs, gs, init, aff, constr, opts)
    outs = {};
    for i = 1:length(opts)
        out = forbes(fs, gs, init, aff, constr, opts{i});
%         disp(out);
        outs{end+1} = out;
    end
end
