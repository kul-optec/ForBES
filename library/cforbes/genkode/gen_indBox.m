function gen_indBox(functionID, upper, lower)
%GEN_INDBOX Generates code for a C function that computes the proximal of
%the indicator function of a box. This function will also compile the
%auto-generated and create a shared library with it.
%
%Syntax:
% gen_indBox(functionID, upper, lower)
%
%Input arguments:
% functionID    : The ID of the function that will be created. The name of
%                 the function is going to be indBox_{functionID} and its
%                 code will be found in indBox_{functionID}.c
% upper, lower  : the upper and lower limits of the box, i.e., the box is
%                 written in the form {x : lower <= x <= upper}. Lower and
%                 upper can either be scalars or vectors. 
%                
%The generated function has the form:
% void indBox_123(double *prox, double *val, const double *x, const int dim);
%
%

functionName = strcat('indBox_',num2str(functionID));
fileName = strcat(functionName,'.c');


fileID = fopen(fileName,'w');       % Open file (Write mode)

% Construct the function
fprintf(fileID,'// Auto-generated C function with ID \n#include "../forbes_helpers.h"\n\n');
fprintf (fileID, 'void %s(\n\tdouble *prox,\n\tdouble *val,\n\tconst double *x,\n\tconst int dim)\n{\n', functionName);
fprintf (fileID, '\tunsigned int i;\n\tfor (i=0; i<dim; i++){\n\t\tprox[i] = MIN(%g, MAX(%g, x[i]));\n\t}\n\t*val = 0.0;\n}', upper, lower);

fclose(fileID);                     % Close the file

% Compile
compile = sprintf('! gcc -c %s.c -o %s.o', ...
    functionName, functionName);
eval(compile);

% Create shared library
shlib = sprintf('! gcc -shared -Wl,-soname,lib%s.so.1 -o lib%s.so.1.0.1  %s.o', ...
    functionName, functionName, functionName);
eval(shlib);

end