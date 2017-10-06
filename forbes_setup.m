% Add ForBES directory to MATLAB's path

forbes_path = fileparts(mfilename('fullpath'));
disp(['Adding ForBES directory to MATLAB path: ', forbes_path]);
addpath(forbes_path);
savepath;

% Compile necessary C source files

utils_path = fullfile(forbes_path, '+forbes', '+utils', filesep);
mex('-outdir', utils_path, [utils_path, 'lbfgs_mex.c'], [utils_path, 'libLBFGS.c']);
mex('-outdir', utils_path, [utils_path, 'RiccatiSolve.c']);

% Success

disp('ForBES was succesfully configured and installed');
disp('Type ''help forbes'' to access the help file');

% Clear variables

clear forbes_path utils_path;
