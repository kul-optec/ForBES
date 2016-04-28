function forbes_test()

forbes_path = fileparts(mfilename('fullpath'));
tests_path = fullfile(forbes_path, 'tests');
addpath(tests_path);

fprintf('%25s', 'CheckGamma... '); test_CheckGamma; fprintf('OK\n');
fprintf('%25s', 'SolveLasso... '); test_SolveLasso; fprintf('OK\n');

rmpath(tests_path);
