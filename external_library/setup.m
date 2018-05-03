% install external libs
lib_path = pwd;

%% install panoc

mex_file_destination = fullfile(lib_path,'forbes_panoc','bin','panoc.mexw64');
if ismac
    disp('Using mex interface Mac os');
    mex_file_location=fullfile(lib_path,'forbes_panoc','bin','panoc_Mac64LLVM.mexw64');
elseif isunix
    disp('Using mex interface Linux');
    mex_file_location=fullfile(lib_path,'forbes_panoc','bin','panoc_Linux64Gcc.mexw64');
elseif ispc
    disp('Using mex interface Windows');
    mex_file_location=fullfile(lib_path,'forbes_panoc','bin','panoc_Windows64VStudio.mexw64');
else
    disp('Platform not supported')
end
copyfile(mex_file_location,mex_file_destination);


%% add paths to path and save
addpath(fullfile(lib_path,'forbes_panoc','bin'));
addpath(fullfile(lib_path,'src_Matlab'));
savepath;% write away path variable

