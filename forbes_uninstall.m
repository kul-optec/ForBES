% Remove ForBES directory from MATLAB's path

forbes_path = fileparts(mfilename('fullpath'));
display(['Removing ForBES directory from MATLAB path: ', forbes_path]);
rmpath(forbes_path);
savepath;

display('ForBES was succesfully removed from MATLAB path');

% Clear variables

clear forbes_path;
