% Remove ForBES directory from MATLAB's path
forbes_path = fileparts(mfilename('fullpath'));
library_path = fullfile(forbes_path, 'library');
private_path = fullfile(forbes_path, 'private');
display(['Removing ForBES directory from MATLAB path: ', forbes_path]);
rmpath(forbes_path);
display(['Removing ForBES library from MATLAB path: ', library_path]);
rmpath(library_path);
savepath;

display('ForBES was succesfully removed from MATLAB path');
