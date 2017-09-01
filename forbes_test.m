function forbes_test()

fprintf('* testing library functions:\n');

% fprintf('%40s', 'sumOp... ');                       tic(); test_sumOp;                        fprintf('OK (%5.2f s)\n', toc());
% fprintf('%40s', 'stackOp... ');                     tic(); test_stackOp;                      fprintf('OK (%5.2f s)\n', toc());

fprintf('%40s', 'separableSum... '); tic();
forbes.tests.test_separableSum;
fprintf('OK (%5.2f s)\n', toc());

fprintf('%40s', 'quadratic... '); tic();
forbes.tests.test_quadratic;
fprintf('OK (%5.2f s)\n', toc());

fprintf('%40s', 'lqrCost... '); tic();
forbes.tests.test_lqrCost;
fprintf('OK (%5.2f s)\n', toc());

fprintf('* testing solver utilities:\n');

fprintf('%40s', 'MakeProblem... '); tic();
forbes.tests.test_MakeProblem;
fprintf('OK (%5.2f s)\n', toc());

fprintf('%40s', 'CheckGamma... '); tic();
forbes.tests.test_CheckGamma;
fprintf('OK (%5.2f s)\n', toc());

fprintf('%40s', 'LineFBE... '); tic();
forbes.tests.test_LineFBE;
fprintf('OK (%5.2f s)\n', toc());

fprintf('%40s', 'SegmentFBE... '); tic();
forbes.tests.test_SegmentFBE;
fprintf('OK (%5.2f s)\n', toc());

% fprintf('%40s', 'FBE inequalities... ');            tic(); forbes.tests.test_inequalities1;                    fprintf('OK (%5.2f s)\n', toc());

fprintf('%40s', 'FBE inequalities... '); tic();
forbes.tests.test_inequalities2;
fprintf('OK (%5.2f s)\n', toc());

% fprintf('* testing composite problems:\n');
%
% fprintf('%36s', 'SolveLasso_small');                tic(); test_SolveLasso_small;             fprintf(' OK (%5.2f s)\n', toc());
% fprintf('%36s', 'SolveLasso_random');               tic(); test_SolveLasso_random;            fprintf(' OK (%5.2f s)\n', toc());
% fprintf('%36s', 'SolveSparseLogReg_small');         tic(); test_SolveSparseLogReg_small ;     fprintf(' OK (%5.2f s)\n', toc());
% fprintf('%36s', 'SolveNuclearNormMC_random');	    tic(); test_SolveNuclearNormMC_random;    fprintf(' OK (%5.2f s)\n', toc());
% % fprintf('%36s', 'SolveRankConstrMC_random');	    tic(); test_SolveRankConstrMC_random;     fprintf(' OK (%5.2f s)\n', toc());
%
% fprintf('* testing separable problems:\n');
%
% fprintf('%36s', 'SolveSVM_random');                 tic(); test_SolveSVM_random;              fprintf(' OK (%5.2f s)\n', toc());
% fprintf('%36s', 'SolveQP_random');                  tic(); test_SolveQP_random;               fprintf(' OK (%5.2f s)\n', toc());
