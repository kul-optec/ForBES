% Add ForBES directory to MATLAB's path

forbes_path = fileparts(mfilename('fullpath'));
disp(['Adding ForBES directory to MATLAB path: ', forbes_path]);
addpath(forbes_path);
savepath;

% Compile necessary C source files

error_msg = 'The C compiler could not succesfully compile ';

LBFGS_dir = fullfile(forbes_path, '+forbes', '+utils', filesep);
LBFGS_src = fullfile(LBFGS_dir, 'lbfgs.c');
if mex('-outdir', LBFGS_dir, LBFGS_src), error([error_msg, LBFGS_src]); end

Riccati_dir = fullfile(forbes_path, '+forbes', '+utils', filesep);
Riccati_src = fullfile(LBFGS_dir, 'RiccatiSolve.c');
if mex('-outdir', Riccati_dir, Riccati_src), error([error_msg, Riccati_src]); end

disp('ForBES was succesfully configured and installed');
disp('Type ''help forbes'' to access the help file');

% Clear variables

clear forbes_path LBFGS_dir LBFGS_src error_msg;
