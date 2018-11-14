function [x,z,y,k,Laug,res_pri,res_dual] = nama_qp(P,q,A,lb,ub,opt,y)

if ~isfield(opt,'maxit')
    opt.maxit = 4000;
end
if ~isfield(opt,'abs_tol')
    opt.abs_tol = 1e-3;
end
if ~isfield(opt,'rel_tol')
    opt.rel_tol = 1e-3;
end

if ~isfield(opt,'memory')
    opt.memory = 10;
end

if ~isfield(opt,'prescale')
    opt.prescale = 1;
end

if ~isfield(opt,'fac')
    opt.fac = 0;
end

if ~isfield(opt,'adaptive')
    opt.adaptive = 0;
end


restart = 0;
pk = [];
qk = [];
cache.cntSkip = 0;
m = size(A,1);


% compute LDL factorization of P
switch opt.prescale
    case 1 %Jacobi scaling
        if opt.fac==1
            LD = ldlchol (P);
            la = sqrt(diag(A*ldlsolve(LD,A')));
            la(la==0)=1;
            e = 1./la;
        else
            L = chol(P);
            la = vecnorm(L\A')';
            la(la==0)=1;
            e = 1./la;
        end
        A = sparse(1:m,1:m,e,m,m)*A;
        lb = e.*lb;ub = e.*ub;
        d=1;c=1;
    case 2 %Ruiz equilibration
        [d,e,c,P,q,A,lb,ub] = RuizEquilibration(P,q,A,lb,ub);
        if opt.fac==1
            LD = ldlchol (P);
        else
            L = chol(P);
        end
    case 3 %Ruiz equilibration to the dual only
        if opt.fac==1
            LD = ldlchol (P);
            Q = A*ldlsolve(LD,A');
            %             c = A*ldlsolve(LD,q);
        else
            L = chol(P);
            Q = L\A'; Q=Q'*Q;
            %             c = A*(L'\(L\q));
        end
        [e,A,lb,ub] = RuizDualEquilibration(Q,A,lb,ub);
        d = 1;c = 1;
    case 0 %Nothing
        if opt.fac==1
            LD = ldlchol (P);
        else
            L = chol(P);
        end
        
        e=1;d = 1; c = 1;
end
if any(isnan(e)) || any(isinf(e))
    keyboard
end

% compute Lipschitz constant and gamma
eigsOpt.issym = 1;
eigsOpt.tol = 1e-3;
if opt.fac == 1
    Lf = eigs(@(y)A*ldlsolve(LD,A'*y),m, 1, 'LM', eigsOpt);
else
    Lf = eigs(@(y)A*(L'\(L\(A'*y))),m, 1, 'LM', eigsOpt);
end

nrmq = norm(q./d,inf);


res_dual = [];
Laug = [];
res_pri = [];
k  = 1;

gam  = 0.95/Lf;
gam = 1.25*gam;

while k < opt.maxit
        Aty = A'*y;
    % x-update x  = -L'\(D\(L\(Aty+q)));
    if opt.fac==1
        x = -ldlsolve(LD,Aty+q);
    else
        x  = -L'\(L\(Aty+q));
    end
    Ax = A*x;
    
    fx = 0.5*(x'*(q-Aty));
    z  = min(max(Ax+y/gam,lb),ub);
    rp = Ax-z;
    Laug(1,k) = fx+y'*rp+0.5*gam*(rp'*rp);
    res_pri(1,k) = norm(rp./e,inf);
    
    if opt.adaptive==1
        y_bar = y + gam*rp;
        Aty_bar   = A'*y_bar;
        if opt.fac==1
            x_bar = -ldlsolve(LD,Aty_bar+q);
        else
            x_bar  = -L'\(L\(Aty_bar+q));
        end
        Ax         = A*x_bar;
        fx_bar     = 0.5*(x_bar'*(q-Aty_bar));
        
        [flag,gam] = check_gamma(x,x_bar,fx,fx_bar,Aty_bar,rp,gam);
        
        if flag == true
            continue;
        end
    end
    
    
    % Compute direction
    [dir, cache] = Direction_lbfgs(opt, k, restart, pk, qk, rp, gam, cache);
    
    % Precompute and store A\*d, inv(P)*A'*d, A*inv(P)*A'*d
    Atd = A'*dir;
    % x-update PiAtd = L'\(D\(L\(Atd)));
    if opt.fac==1
        PiAtd = ldlsolve(LD,Atd);
    else
        PiAtd = L'\(L\(Atd));
    end
    APiAtd = A*PiAtd;
    % Try unit stepsize first
    yt    = y+dir;
    Atyt  = Aty+Atd;
    xt    = x-PiAtd;
    Axt   = Ax-APiAtd;
    zt    = min(max(Axt+yt/gam,lb),ub);
    fxt   = 0.5*(xt'*(q-Atyt));
    rpt   = Axt-zt;
    Laugt = fxt+yt'*rpt+0.5*gam*(rpt'*rpt);
    tau   = 1;
    
    if Laugt < Laug(1,k)
        Atrp = A'*rp;
        if res_pri(1,k) <= opt.abs_tol+opt.rel_tol*max(norm(Ax./e,inf),norm(z./e,inf))
            nrmAty = Aty+gam*Atrp;
            nrmPx  = norm((nrmAty+q)./d,inf);
            nrmAty = norm(nrmAty./d,inf);
            res_dual(1,k) = gam*norm(Atrp./d,inf)/c;
            if res_dual(1,k) <= opt.abs_tol+opt.rel_tol*max([nrmPx,nrmAty,nrmq])
                x = d.*x; z = z./e; y = (e.*(y+gam*rp))/c;
                break
            end
        end
        
        % PiAtrp = L'\(D\(L\(Atrp)));
        if opt.fac==1
            PiAtrp = ldlsolve(LD,Atrp);
        else
            PiAtrp = L'\(L\(Atrp));
        end
        APiAtrp = A*PiAtrp;
        tau_ctr = 0;
        
        while Laugt <= Laug(1,k)
            tau   = 0.5*tau;
            if tau <= 1e-3
                break
            end
            %yt    = (1 - tau)*y_bar + tau*(y+dir);
            yt    = y+tau*dir+(gam*(1-tau))*rp;
            Atyt  = Aty+tau*Atd+(gam*(1-tau))*Atrp;
            xt    = x-tau*PiAtd-(gam*(1-tau))*PiAtrp;
            Axt   = Ax-tau*APiAtd-(gam*(1-tau))*APiAtrp;
            zt    = min(max(Axt+yt/gam,lb),ub);
            fxt   = 0.5*(xt'*(q-Atyt));
            rpt   = Axt-zt;
            Laugt = fxt+yt'*rpt+0.5*gam*(rpt'*rpt);
            tau_ctr = tau_ctr + 1;
        end
        
    end
    
    % store pairs for LBFGS
    pk   = yt - y;
    %     qk = (Ax-min(max(Ax+y/10,lb),ub))-(Axt-min(max(Axt+yt/10,lb),ub));
    
    
    while 1
        zt    = min(max(Axt+yt/gam,lb),ub);
        rpt   = Axt-zt;
        qk = rp - rpt;
        
        if opt.adaptive
            [flag,gam] = check_gamma(xt,x,fxt,fx,Aty,rpt,gam);
            if ~flag
                break
            end
        else
            break
        end
    end
    
    y = yt+gam*rpt;
    if mod(k,0) == 0
        Aty = A'*y;
    else
        Atrpt = A'*rpt;
        Aty = Atyt+gam*Atrpt;
        if norm(rpt./e,inf) <= opt.abs_tol+opt.rel_tol*max(norm(Axt./e,inf),norm(zt./e,inf))
            res_dual(1,k) = gam*norm(Atrpt./d,inf)/c;
            if res_dual(1,k) <= opt.abs_tol+opt.rel_tol*max([norm((Atyt+q)./d,inf),norm(Aty./d,inf),nrmq])
                x = d.*xt; z = zt./e;y= (e.*y)/c;
                break
            end
        end
    end
    % x = -L'\(D\(L\(Aty+q)));
    if opt.fac == 1
        x = -ldlsolve(LD,Aty+q);
    else
        x = -L'\(L\(Aty+q));
    end
    
    Ax = A*x;
    if norm((Ax-zt)./e,inf)<= opt.abs_tol+opt.rel_tol*max(norm(Ax./e,inf),norm(zt./e,inf))
        x = d.*x; z = zt./e; y = (e.*y)/c;
        
        break
    end
    z  = min(max(Ax+y/gam,lb),ub);
    rp = Ax-z;
    fx = 0.5*(x'*(q-Aty));
    k = k+1;
    Laug(1,k+1)    = fx+y'*rp+0.5*gam*(rp'*rp);
    res_pri(1,k+1) = norm(rp./e,inf);
    
    if opt.adaptive==1
        gam = 1.25*gam;   
    end
end
end
% L-BFGS

function [dir, cache] = Direction_lbfgs(opt, it, restart, sk, yk, rp, gam, cache)

sk = full(sk(:));
yk = full(yk(:));

if it == 1 || restart
    dir = gam*rp; % use steepest descent direction initially
    cache.LBFGS_col = 0; % last column of Sk, Yk that was filled in
    cache.LBFGS_mem = 0; % current memory of the method
else
    YSk = yk'*sk;
    if YSk > 0.0
        cache.LBFGS_col = 1+mod(cache.LBFGS_col, opt.memory);
        cache.LBFGS_mem = min(cache.LBFGS_mem+1, opt.memory);
        cache.S(:,cache.LBFGS_col) = sk;
        cache.Y(:,cache.LBFGS_col) = yk;
        cache.YS(cache.LBFGS_col) = YSk;
    else
        cache.cntSkip = cache.cntSkip+1;
    end
    if cache.LBFGS_mem > 0
        H = cache.YS(cache.LBFGS_col)/...
            (cache.Y(:,cache.LBFGS_col)'*cache.Y(:,cache.LBFGS_col));
        dir = lbfgs(cache.S, cache.Y, cache.YS, H, ...
            rp, int32(cache.LBFGS_col), int32(cache.LBFGS_mem));
    else
        dir = gam*rp;
    end
end

end



%  new_gam = adaptive(gam, res_pri, res_dual, 2, 100);
%  gam = min(gam, new_gam);
% function new_gam = adaptive(gam, pri_res, dual_res, mu, tau)
% if norm(pri_res) > mu * norm(dual_res)
%     new_gam = tau * gam;
% elseif norm(dual_res) > mu * norm(pri_res)
%     new_gam = 1/tau * gam;
% else
%     new_gam = gam;
% end
% end
%
