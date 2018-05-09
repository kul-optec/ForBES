function [ out ] = solve_panoc( f, g ,opt , intial_solution )
    %SOLVE_PANOC Solve problem with panoc lib
    
    if(~isfield(g,'name'))
        g.name='costum';
    end
    
    % parse the options into the panoc format
    problem.dimension = length(intial_solution);
    
    if(strcmp(g.name,'box'))
        problem.constraint_type = 'box';
        
        problem.upper_bound=g.upper_bound;
        problem.lower_bound=g.lower_bound;
    else % unknown function, use the matlab version provided by the user
        disp(['WARNING: No C implementation found of the constraint ' ... 
                'function using the matlab implementation,'...
                ' this may impact performance.']);
        problem.constraint_type = 'costum';
        
        % convert the function into the right format and set gamma=0
        g_function = g.makeprox();
        constraint = @(x,gamma) g_function(x,gamma);
        problem.constraint = constraint;
    end

    solver_params.tolerance = opt.tol;
    solver_params.buffer_size = 20;
    solver_params.max_iterations = opt.maxit;
    
    % Solve the problem with panoc,
    % copy over by value the initial position.
    solution=zeros(problem.dimension,1);
    for i=1:problem.dimension
        solution(i)=intial_solution(i);
    end
    
    panoc('init',problem,solver_params);
    number_of_iterations = panoc('solve',solution,f);
    panoc('cleanup');
    
    out.x=solution;
    out.iterations=number_of_iterations;
end