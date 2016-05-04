function [ S ] = system_masses(N, options )
%masses_model function is used to find the dynamic equations of N masses
%system. The no of masses is specified by N and options specify the constraints
%on the system.
% The syntax of the function is
% S=masses_model(N,options)
% system with N-1 inputs .
default_options = struct('M', 1*ones(N,1), 'b', 0.1*ones(N+1,1),...
    'k',1*ones(N+1,1),'xmin',-5*ones(2*N,1), 'xmax', 5*ones(2*N,1), 'umin', ...
    -5*ones(N-1,1),'umax',5*ones(N-1,1), 'Ts', 0.1);
ops = default_options;
if nargin==2, % User-provided options
    ops = options;
    flds = fieldnames(default_options);
    for i=1:numel(flds),
        if ~isfield(options, flds(i))
            ops.(flds{i})=default_options.(flds{i});
        end
    end
end

M=ops.M;        % M(i)=mass of body #i
b=ops.b;        % b(i)=viscous friction of body #i
k=ops.k;        % k(i)=spring of body #i


% Define full A,B model
nx=2*N;
nu=N-1;
S.nx=nx;
S.nu=nu;


Ag=zeros(nx,nx);
Bg=zeros(nx,nu);
for i=1:N,
    h=2*i-1;
    Ag(h,h+1)=1;  % velocity
    Ag(h+1,h+1)=-(b(i)+b(i+1))/M(i);  % friction
    Ag(h+1,h)=-(k(i)+k(i+1))/M(i);    % self-springs
    if i>1,
        Ag(h+1,h-2)=k(i)/M(i);
        Ag(h+1,h-1)=b(i)/M(i);
    end
    if i<N,
        Ag(h+1,h+2)=k(i+1)/M(i);
        Ag(h+1,h+1)=b(i+1)/M(i);
    end
    if i>1
        Bg(h+1,i-1)=-1/M(i);
    end
    if i<N
        Bg(h+1,i)=1/M(i);
    end
end
%Bg
%S.A=eye(2*N)+ops.Ts*Ag;
%S.B=ops.Ts*Bg;

sysgc=ss(Ag,Bg,eye(nx),zeros(nx,nu));

% Discretize and normalise...
sysgd=c2d(sysgc, ops.Ts);
tol=1e-6;
S.A=(1-(abs(sysgd.a)<tol)).*sysgd.a;
S.B=(1-(abs(sysgd.b)<tol)).*sysgd.b;
S.F=[eye(S.nx);-eye(S.nx);zeros(2*S.nu,S.nx)];
S.G=[zeros(2*S.nx,S.nu);eye(S.nu);-eye(S.nu)];
S.g=[ops.xmax;-ops.xmin;ops.umax;-ops.umin];

end
