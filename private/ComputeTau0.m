function tau0 = ComputeTau0(prob, opt, it, hasGammaChanged, sk, yk, dir)

switch opt.methodID
    case {2, 3, 8, 9, 10} % (limited-memory) quasi-Newton methods
        tau0 = 1.0;
    case 7 % Barzilai-Borwein
        if it == 1 || hasGammaChanged
            tau0 = 1.0/norm(dir, inf);
        else
            tau0 = (sk'*sk)/(sk'*yk);
        end
end

end
