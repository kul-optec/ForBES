% Only use this demo after installing the panoc library. This is only
% supported on Windows right now. Run ./external_library/setup.m to
% install panoc.

clear all;
%%
f=@rosen;
g=indBox(-4,4);
x0=[0;0];
aff=0; % not supported with panoc
opts.solver='panoc';
opts.tol=1e-12;
opts.maxit=200;

% g.name='blah'; % uncomment this line if you want to use the Matlab function constraint
tic
out = forbes(f, g, x0, aff, [], opts) % this problem should take 20 iterations
toc

% double check the results
if(out.iterations~=20)
    disp(['Error: not solved in 20 iterations as expected but in ' ...
        num2str(out.iterations) 'iterations'] );
end

theoretical_solution = [1;1];
if(norm(out.x-theoretical_solution)>opts.tol)
    disp(['Error: solution not [1;1] as expected but[' ...
        num2str(out.x(1)) ';' num2str(out.x(2)) ']' ]);
end