function [diagS,c,Pbar,qbar,E,D,Abar,lbar,ubar] = ruiz_equilibration(P,q,A,l,u)
eps_equi = 1e-1;
c = 1;
[m,n] = size(A);
mn = m+n; 
diagS = ones(mn,1);
delta = zeros(mn,1);
Pbar = P; Abar = A; qbar = q;lbar = l; ubar = u;
while norm(delta-1,inf)>eps_equi
%     colPnorm = vecnorm(Pbar, 'inf');
%     colAnorm = vecnorm(Abar, 'inf');
%     colA_tnorm = vecnorm(Abar', 'inf');
%     vec_norms = [max(colPnorm, colAnorm)'; colA_tnorm'];
%     delta = 1./sqrt(vec_norms);
    delta = 1./sqrt(vecnorm([Pbar Abar';Abar sparse(m,m)],'inf')');
    d = delta(1:n,1);
    e = delta(n+1:mn,1);
    D = sparse(1:n,1:n,d,n,n);
    E = sparse(1:m,1:m,e,m,m);
    Pbar = D*Pbar*D;
    qbar = d.*qbar;
    Abar = E*Abar*D;
    lbar = e.*lbar;
    ubar = e.*ubar;
    gam = 1/max(mean(vecnorm(Pbar,'inf')),norm(qbar,'inf'));
    Pbar = gam*Pbar;
    qbar = gam*qbar;
    diagS = delta.*diagS;
    c = gam*c;
end
