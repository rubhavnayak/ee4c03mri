clc; clear; close all;

slice = 8;
[badCh1, badCh2, badCh3, goodCh1, goodCh2, goodCh3] = loadSlice(slice);
fprintf('\nSlice No: %d\n', slice)

% Reconstruct images from individual channels
img1 = kSpaceToImage(goodCh1, zeros(size(goodCh1)), zeros(size(goodCh1)));
img2 = kSpaceToImage(zeros(size(goodCh2)), goodCh2, zeros(size(goodCh2)));
img3 = kSpaceToImage(zeros(size(goodCh3)), zeros(size(goodCh3)), goodCh3);

figure;
subplot(1, 3, 1); imagesc(img1); colormap gray; title('Channel 1');
subplot(1, 3, 2); imagesc(img2); colormap gray; title('Channel 2');
subplot(1, 3, 3); imagesc(img3); colormap gray; title('Channel 3');

% Calculate variance in top 10% of rows for each channel as noise estimate
noise_var_ch1 = var(goodCh1(1:50, :), 0, 'all');
noise_var_ch2 = var(goodCh2(1:50, :), 0, 'all');
noise_var_ch3 = var(goodCh3(1:50, :), 0, 'all');

fprintf('Noise Variance - Channel 1: %.4f\n', noise_var_ch1);
fprintf('Noise Variance - Channel 2: %.4f\n', noise_var_ch2);
fprintf('Noise Variance - Channel 3: %.4f\n', noise_var_ch3);

% Reconstruct images from each channel
img1 = kSpaceToImage(goodCh1, zeros(size(goodCh1)), zeros(size(goodCh1)));
img2 = kSpaceToImage(zeros(size(goodCh2)), goodCh2, zeros(size(goodCh2)));
img3 = kSpaceToImage(zeros(size(goodCh3)), zeros(size(goodCh3)), goodCh3);
imgK = kSpaceToImage(goodCh1, goodCh2, goodCh3);

img1 = double(img1);
img2 = double(img2);
img3 = double(img3);

% Perform root sum of squares fusion
fused_img = sqrt(abs(img1).^2 + abs(img2).^2 + abs(img3).^2);

figure;
subplot(2, 2, 1); imagesc(img1); colormap gray; title('Channel 1');
subplot(2, 2, 2); imagesc(img2); colormap gray; title('Channel 2');
subplot(2, 2, 3); imagesc(img3); colormap gray; title('Channel 3');
subplot(2, 2, 4); imagesc(imgK); colormap gray; title('Using KSpacetoImage()');
