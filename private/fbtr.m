function out = fbtr(prob, opt, lsopt)

% initialize output stuff

if opt.report
    residual = zeros(1, opt.maxit);
    objective = zeros(1, opt.maxit);
    ts = zeros(1, opt.maxit);
    % initialize operations counter
    ops = FBOperations();
else
    ops = [];
end

% get Lipschitz constant & adaptiveness

[Lf, adaptive] = prob.Get_Lipschitz(opt);

% initialize gamma

gam = (1-opt.beta)/Lf;

% display header

if opt.display >= 2
    fprintf('\n%s', opt.name);
    fprintf('\n%6s%11s%11s%11s%11s%11s%11s\n', 'iter', 'gamma', 'optim.', 'object.', '||dir||', 'slope', 'tau');
end

cache_dir.cntSkip = 0;
restart = 0;
delta_k = opt.delta_k;
Hk = opt.Hk;
Bk = opt.Bk;
cache_initial = FBCache(prob, prob.x0, gam, ops);

t0 = tic();

for it = 1:opt.maxit

    % store initial cache

    if it == 1
        cache_0 = cache_initial;
        cache_previous = cache_initial;
    end

    % trace stuff

    prox_point    = cache_previous.Get_ProxGradStep();
    cache_current = FBCache(prob, prox_point, gam, ops);
    
    if opt.report
        objective(1, it) = cache_current.Get_FBE();
        ts(1, it) = toc(t0);
        residual(1, it) = norm(cache_current.Get_FPR(), 'inf')/cache_current.Get_Gamma();
    end
    if opt.toRecord
        record(:, it) = opt.record(prob, it, cache_0, cache_current);
    end

    % check for termination

    if ~restart
        if ~opt.customTerm
            if cache_current.Check_StoppingCriterion(opt.tol)
                msgTerm = 'reached optimum (up to tolerance)';
                flagTerm = 0;
                break;
            end
        else
            %flagStop = opt.term(prob, it, cache_0, cache_current);
            flagStop = opt.term(prob, it, cache_0, cache_current,opt.hack,opt.H);
            if (adaptive == 0 || it > 1) && flagStop
                msgTerm = 'reached optimum (custom criterion)';
                flagTerm = 0;
                break;
            end
        end
    end
  
     if it == 1 || restart 
        sk = [];
        yk = [];
     end
   
    dir             = opt.subProb(Hk, Bk, cache_current.Get_GradFBE(), delta_k);
    cache_current   =  FBCache(prob, cache_current.x + dir, gam, ops);
    ared            =  cache_previous.Get_FBE()-cache_current.Get_FBE();
    pred            = -(cache_current.Get_GradFBE()'*dir + 0.5*dir'*Bk*dir);
    rkh             =  ared/pred;
    p               =  0;
    max_inner       =  100;
   while (rkh < 0.25 && p < max_inner)
       p               =  p + 1;
       delta_k         =  0.5*delta_k;
       dir             =  opt.subProb(Hk, Bk, cache_current.Get_GradFBE(), delta_k);
       cache_current   =  FBCache(prob, cache_current.x + dir, gam, ops);
       ared            =  cache_previous.Get_FBE()-cache_current.Get_FBE();
       pred            = -(cache_current.Get_GradFBE()'*dir + 0.5*dir'*Bk*dir);
       rkh             =  ared/pred;
   end
    pk  = p;
    %fprintf('Number of inner while loop iterations %5d\n',pk);   
   
    if (rkh >= 0.75 && norm(dir,2)>((1-eps)*delta_k))     
        delta_k     = 1.5*delta_k;
        
    end
   
    %% store pair (s, y) to compute next direction
    
    sk  = cache_current.Get_Point()    - cache_previous.Get_Point();
    yk  = cache_current.Get_GradFBE()  - cache_previous.Get_GradFBE();

    %cache_current  = cache_current.Get_ProxGradStep();       
    Sflag          = 1;
    Bk             = bfgs(Bk,sk,yk,Sflag);
    Sflag          = 2; 
    Hk             = bfgs(Hk,sk,yk,Sflag);
    
    % prepare next iteration, store current solution
    cache_previous = cache_current;
 
    % display stuff

    if opt.display == 1
        Util_PrintProgress(it);
    elseif (opt.display == 2 && mod(it,10) == 0) || opt.display >= 3
        fprintf('%6d %7.4e %7.4e %7.4e %7.4e \n', it, gam, residual(1,it), objective(1,it), norm(dir));
    end

end

time = toc(t0);

if it == opt.maxit
    flagTerm = 1;
    msgTerm = 'exceeded maximum iterations';
end

if opt.display == 1
    Util_PrintProgress(it, flagTerm);
elseif opt.display >= 2
    fprintf('%6d %7.4e %7.4e %7.4e\n', it, gam, residual(1,it), objective(1,it));
end

% pack up results

out.name = opt.name;
out.message = msgTerm;
out.flag = flagTerm;
if it == opt.maxit
    out.x = cache_current.Get_Point();
else
    out.x = cache_current.Get_ProxGradStep();
end
out.iterations = it;
out.operations = ops;
if opt.report
    out.residual = residual(1, 1:it);
    out.objective = objective(1, 1:it);
    out.ts = ts(1, 1:it);
end
if opt.toRecord, out.record = record; end
out.gam = gam;
out.time = time;
out.adaptive = adaptive;
end
