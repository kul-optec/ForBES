function forbes_test()

% build the paths to add
forbes_path = fileparts(mfilename('fullpath'));
tests_path = fullfile(forbes_path, 'tests');
addpath(tests_path);

fprintf('%40s', 'CheckGamma... '); test_CheckGamma; fprintf('OK\n');
fprintf('%40s', 'SolveLasso... '); test_SolveLasso; fprintf('OK\n');
fprintf('%40s', 'SolveSparseLogReg... '); test_SolveSparseLogReg; fprintf('OK\n');
fprintf('%40s', 'SVM... '); test_SolveSVM; fprintf('OK\n');
fprintf('%40s', 'NuclearNormMC... '); test_SolveNuclearNormMC; fprintf('OK\n');
fprintf('%40s', 'RankConstrMC... '); test_SolveRankConstrMC; fprintf('OK\n');

% rebuild the paths to remove (test scripts may have cleared)
forbes_path = fileparts(mfilename('fullpath'));
tests_path = fullfile(forbes_path, 'tests');
rmpath(tests_path);
