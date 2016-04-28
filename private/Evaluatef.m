function [v, cnt] = Evaluatef(prob, x)
    %      Q,C1,C2,f2, g
    cnt = [0, 0, 0, 0, 0];
    f1x = 0; f2x = 0;
    if prob.istheref1
        if prob.isthereC1
            if prob.isC1fun, C1x = prob.C1(x);
            else C1x = prob.C1*x; end
            res1x = C1x - prob.d1;
            if prob.isQfun, Qres1x = prob.Q(res1x);
            else Qres1x = prob.Q*res1x; end
            cnt(2) = cnt(2)+1;
        else
            res1x = x - prob.d1;
            if prob.isQfun, Qres1x = prob.Q(res1x);
            else Qres1x = prob.Q*res1x; end
        end
        cnt(1) = cnt(1)+1;
        f1x = 0.5*(res1x'*Qres1x) + prob.q'*res1x;
    end
    if prob.istheref2
        if prob.isthereC2
            if prob.isC2fun, C2x = prob.C2(x);
            else C2x = prob.C2*x; end
            res2x = C2x - prob.d2;
            f2x = prob.callf2(res2x);
            cnt(3) = cnt(3)+1;
        else
            res2x = x - prob.d2;
            f2x = prob.callf2(res2x);
        end
        cnt(4) = cnt(4)+1;
    end
    if prob.istherelin
        v = f1x + f2x + prob.l'*x;
    else
        v = f1x + f2x;
    end
end