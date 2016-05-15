function forbes_test()

% build the paths to add
forbes_path = fileparts(mfilename('fullpath'));
tests_path = fullfile(forbes_path, 'tests');
addpath(tests_path);

fprintf('%40s', 'CheckGamma... ');                  test_CheckGamma;                   fprintf('OK\n');
fprintf('%40s', 'DirFBE... ');                      test_DirFBE;                       fprintf('OK\n');
fprintf('%40s', 'SolveLasso_small... ');            test_SolveLasso_small;             fprintf('OK\n');
fprintf('%40s', 'SolveLasso_random... ');           test_SolveLasso_random;            fprintf('OK\n');
fprintf('%40s', 'SolveSparseLogReg_small... ');     test_SolveSparseLogReg_small ;     fprintf('OK\n');
fprintf('%40s', 'SolveSVM_random... ');             test_SolveSVM_random;              fprintf('OK\n');
fprintf('%40s', 'SolveNuclearNormMC_random... ');	test_SolveNuclearNormMC_random;    fprintf('OK\n');
fprintf('%40s', 'SolveRankConstrMC_random... ');	test_SolveRankConstrMC_random;     fprintf('OK\n');

% rebuild the paths to remove (test scripts may have cleared)
forbes_path = fileparts(mfilename('fullpath'));
tests_path = fullfile(forbes_path, 'tests');
rmpath(tests_path);
