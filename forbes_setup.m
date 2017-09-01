% Add ForBES directory to MATLAB's path
forbes_path = fileparts(mfilename('fullpath'));
disp(['Adding ForBES directory to MATLAB path: ', forbes_path]);
addpath(forbes_path);
savepath;

% Compile necessary C source files
LBFGS_path = fullfile(forbes_path, 'private', 'lbfgs.c');
error_msg = 'The C compiler could not succesfully compile ';
if mex('-outdir', private_path, LBFGS_path), error([error_msg, LBFGS_path]); end
disp('ForBES was succesfully configured and installed');
disp('Type ''help forbes'' to access the help file');

% Clear variables
clear forbes_path LBFGS_path error_msg;
