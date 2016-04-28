% solve a matrix decomposition problem for background/foreground extraction

close all;
clear;

rng(0);

% customize paths and other things

basedir = './ShoppingMall/';
basename = 'ShoppingMall';
extension = 'bmp';
idx_start = 1060;
idx_end = 1109;

% load the frames
IMGS = []; % stack the frames into a vector
for idx = idx_start:idx_end
    filename = [basedir, basename, num2str(idx), '.', extension];
    IMG = imread(filename);
    if length(size(IMG)) == 3
        IMG = rgb2gray(IMG);
    end
    H = size(IMG, 1);
    W = size(IMG, 2);
    M = H*W; % number of pixels per frame
    IMGS = [IMGS; double(IMG(:))/255];
end
N = idx_end-idx_start+1; % number of frames
IMGS_MAT = reshape(IMGS, M, N);

% set parameters

r = 1;
lam = 1e-2; % l0 norm coefficient

% run forbes

S = [speye(M*N), speye(M*N)]; % linear operator that sums two vectors
f = quadLoss(1, IMGS);
g = separableSum({indRankBall(M, N, r), l0Norm(lam)}, [M*N, M*N]);

opt.display = 2;
opt.Lf = 2;
opt.tol = 1e-8;
opt.record = @(prob, it, gam, c0, c, ops) [ops.proxg];

opt.method = 'fbs'; opt.variant = 'basic';
out_fbs = forbes(f, g, [IMGS; zeros(M*N, 1)], S, {}, opt);
x_back = out_fbs.x(1:M*N);
x_fore = out_fbs.x((M*N+1):end);
fprintf('FBS\n');
fprintf('Time             : %.2f sec\n', out_fbs.ts(end));
fprintf('Iterations       : %d\n', out_fbs.iterations);
fprintf('SVDs             : %d\n', out_fbs.operations.proxg);
fprintf('Residual (rel)   : %7.4e\n', norm(S*out_fbs.x-IMGS)/norm(IMGS));
fprintf('Rank     (X1)    : %d\n', rank(reshape(x_back,M,N)));
fprintf('L1 norm  (X2)    : %7.4e\n', norm(x_fore,1));
fprintf('NNZ      (X2)    : %d\n', nnz(x_fore));

opt.method = 'lbfgs-fpr-old'; opt.variant=''; opt.memory = 10;
out_lbfgs_fpr = forbes(f, g, [IMGS; zeros(M*N, 1)], S, {}, opt);
out = out_lbfgs_fpr;
x_back = out.x(1:M*N);
x_fore = out.x((M*N+1):end);
fprintf('L-BFGS\n');
fprintf('Time             : %.2f sec\n', out.ts(end));
fprintf('Iterations       : %d\n', out.iterations);
fprintf('SVDs             : %d\n', out.operations.proxg);
fprintf('Residual (rel)   : %7.4e\n', norm(S*out.x-IMGS)/norm(IMGS));
fprintf('Rank     (X1)    : %d\n', rank(reshape(x_back,M,N)));
fprintf('L1 norm  (X2)    : %7.4e\n', norm(x_fore,1));
fprintf('NNZ      (X2)    : %d\n', nnz(x_fore));

% save result images

for idx = 0:N-1
    filename_bw = [basedir, 'out/', 'ORIG_', basename, num2str(idx_start+idx), '.png'];
    filename_back = [basedir, 'out/', 'BACK_', basename, num2str(idx_start+idx), '.png'];
    filename_fore = [basedir, 'out/', 'FORE_', basename, num2str(idx_start+idx), '.png'];
    x_bw_curr = reshape(IMGS(idx*M+1:(idx+1)*M), H, W);
    x_back_curr = reshape(x_back(idx*M+1:(idx+1)*M), H, W);
    x_fore_curr = reshape(x_fore(idx*M+1:(idx+1)*M), H, W);
%     figure(1); imshow(x_bw_curr);
%     figure(2); imshow(x_back_curr);
%     figure(3); imshow(x_fore_curr);
    imwrite(x_bw_curr, filename_bw);
    imwrite(x_back_curr, filename_back);
    imwrite(x_fore_curr, filename_fore);
end
