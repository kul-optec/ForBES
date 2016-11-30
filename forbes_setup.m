% Add ForBES directory to MATLAB's path
forbes_path = fileparts(mfilename('fullpath'));
library_path = fullfile(forbes_path, 'library');
private_path = fullfile(forbes_path, 'private');
display(['Adding ForBES directory to MATLAB path: ', forbes_path]);
addpath(forbes_path);
display(['Adding ForBES library to MATLAB path: ', library_path]);
addpath(library_path);
savepath;

% Clear variables
clear forbes_path, library_path, private_path;

% Compile necessary C source files
LBFGS_path = fullfile(forbes_path, 'private', 'lbfgs.c');
Riccati_path = fullfile(forbes_path, 'library', 'RiccatiSolve.c');
error_msg = 'The C compiler could not succesfully compile ';
if mex('-outdir', private_path, LBFGS_path), error([error_msg, LBFGS_path]); end
if mex('-outdir', library_path, Riccati_path), error([error_msg, Riccati_path]); end
display('ForBES was succesfully configured and installed');
display('Type ''help forbes'' to access the help file');

% Clear variables
clear LBFGS_path, Riccati_path, error_msg;
