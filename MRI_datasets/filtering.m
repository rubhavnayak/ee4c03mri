%% CLEAR OUTPUTS AND WORKSPACE, INITIALIZE FUNCTIONS
clc; clear; close all;

%Median Filter
%Arguments: Three channels of K-space images, vector of outlier columns
function [ch1_filtered, ch2_filtered, ch3_filtered, filtered_eye] = medianFilter(ch1, ch2, ch3, good_eye, outlier_columns)
    
    bad_eye = kSpaceToImage(ch1, ch2, ch3);
    mse_good_bad = mean((good_eye(:) - bad_eye(:)).^2);
    fprintf('MSE between Reference and Original Image: %.4f\n', mse_good_bad);
   
    p = 2;
    filtered_channels = {ch1, ch2, ch3};

    for idx = 1:3
        channel = filtered_channels{idx};
        
        for col_idx = outlier_columns
            window_start = max(col_idx - p, 1);
            window_end = min(col_idx + p, size(channel, 2));
            window_data = channel(:, window_start:window_end);
            
            % Apply median filter across columns in the window
            estimated_column = median(window_data, 2);
            
            % Replace corrupted column with median-filtered data
            channel(:, col_idx) = estimated_column;
        end
        
        filtered_channels{idx} = channel;
    end

    % Output corrected channels
    ch1_filtered = filtered_channels{1};
    ch2_filtered = filtered_channels{2};
    ch3_filtered = filtered_channels{3};


    filtered_eye = kSpaceToImage(ch1_filtered, ch2_filtered, ch3_filtered);
    mse_good_filt = mean((good_eye(:) - filtered_eye(:)).^2);
    fprintf('MSE between Reference and Filtered Image: %.4f\n', mse_good_filt);

    mse_diff = mse_good_bad - mse_good_filt;
    fprintf('MSE Diff: %.4f\n', mse_diff);
end

%Wiener Filter Implemented From Lecture
%Arguments: Three channels of K-space images, vector of outlier columns, window size
function [ch1_filtered, ch2_filtered, ch3_filtered] = WienerFilter(ch1, ch2, ch3, outlier_columns, window_size)

    filtered_channels = {ch1, ch2, ch3};

    for idx = 1:3
        channel = filtered_channels{idx};
   
        noise_variance = var(channel(:));

        for col_idx = outlier_columns
            % Adjacent Columns
            window_start = max(col_idx - window_size, 1);
            window_end = min(col_idx + window_size, size(channel, 2));
            window_data = channel(:, window_start:window_end);
            
            
            local_mean = mean(window_data, 2);
            local_variance = var(window_data, 0, 2);
            
            % Apply Wiener filter formula
            estimated_column = local_mean + ...
                (max(local_variance - noise_variance, 0) ./ max(local_variance, noise_variance)) .* ...
                (window_data(:, ceil(end/2)) - local_mean);
            
            channel(:, col_idx) = estimated_column;
        end
        
        filtered_channels{idx} = channel;
    end

    ch1_filtered = filtered_channels{1};
    ch2_filtered = filtered_channels{2};
    ch3_filtered = filtered_channels{3};

end

%Wiener #2, using matlab wiener2()
%Arguments: Three channels of K-space images, vector of outlier columns, window size
function [filteredCh1, filteredCh2, filteredCh3] = WienerFilter2(kSpaceCh1, kSpaceCh2, kSpaceCh3, outlierCols, window)
  
    filteredCh1 = kSpaceCh1;
    filteredCh2 = kSpaceCh2;
    filteredCh3 = kSpaceCh3;
    filterWindow = [window, 1];

    for i = 1:length(outlierCols)
        colIdx = outlierCols(i);

        % Apply Wiener filter to the real and imaginary parts separately
        filteredCh1(:, colIdx) = wiener2(real(kSpaceCh1(:, colIdx)), filterWindow) + ...
                                  1i * wiener2(imag(kSpaceCh1(:, colIdx)), filterWindow);
        filteredCh2(:, colIdx) = wiener2(real(kSpaceCh2(:, colIdx)), filterWindow) + ...
                                  1i * wiener2(imag(kSpaceCh2(:, colIdx)), filterWindow);
        filteredCh3(:, colIdx) = wiener2(real(kSpaceCh3(:, colIdx)), filterWindow) + ...
                                  1i * wiener2(imag(kSpaceCh3(:, colIdx)), filterWindow);
    end
end

% Adaptive Wiener
%Arguments: Three channels of K-space images, vector of outlier columns, window size
function [filteredCh1, filteredCh2, filteredCh3] = WienerFilter2Adaptive(kSpaceCh1, kSpaceCh2, kSpaceCh3, outlierCols, window)
    
    filteredCh1 = kSpaceCh1;
    filteredCh2 = kSpaceCh2;
    filteredCh3 = kSpaceCh3;

    for i = 1:length(outlierCols)
        colIdx = outlierCols(i);
        
        % Calculate the local variance for each channel in the column
        localVariance = var([real(kSpaceCh1(:, colIdx)); real(kSpaceCh2(:, colIdx)); real(kSpaceCh3(:, colIdx))]);
        %disp(localVariance);
        
        % Define adaptive window size
        if localVariance < 10 %0.05  
            filterWindow = [window, 1];
        elseif localVariance < 100 %0.2
            filterWindow = [window + 10, 1];
        else
            filterWindow = [window + 20, 1];
        end
        
        % Apply Wiener filter 
        filteredCh1(:, colIdx) = wiener2(real(kSpaceCh1(:, colIdx)), filterWindow) + ...
                                  1i * wiener2(imag(kSpaceCh1(:, colIdx)), filterWindow);
        filteredCh2(:, colIdx) = wiener2(real(kSpaceCh2(:, colIdx)), filterWindow) + ...
                                  1i * wiener2(imag(kSpaceCh2(:, colIdx)), filterWindow);
        filteredCh3(:, colIdx) = wiener2(real(kSpaceCh3(:, colIdx)), filterWindow) + ...
                                  1i * wiener2(imag(kSpaceCh3(:, colIdx)), filterWindow);
    end
end

% Best Wiener Filter (Automatically computes best window size, uses WienerFilter())
%Arguments: Three channels of K-space images, reconstructed good image, vector of outlier columns
function [ch1_filtered, ch2_filtered, ch3_filtered, filtered_eye] = BestWienerFilter(badCh1, badCh2, badCh3, good_eye, outlier_columns)
   
    %Get reference MSE
    bad_eye = kSpaceToImage(badCh1, badCh2, badCh3);
    mse_good_bad = mean((good_eye(:) - bad_eye(:)).^2);
    fprintf('MSE between Reference and Original Image: %.4f\n', mse_good_bad);

    % Initialize variables
    best_mse_diff = -Inf;
    best_window_size = 0;

    % Loop through window sizes
    for i = 2:250
        
        window_size = 1 + (i - 1) * 2;
        
        % Apply Wiener filter with the current window size
        %[ch1, ch2, ch3] = WienerFilter(badCh1, badCh2, badCh3, outlier_columns, window_size);
        %[ch1, ch2, ch3] = WienerFilter2(badCh1, badCh2, badCh3, outlier_columns, window_size);
        [ch1, ch2, ch3] = WienerFilter2Adaptive(badCh1, badCh2, badCh3, outlier_columns, window_size);
 
        filtered_eye = kSpaceToImage(ch1, ch2, ch3);
        
        mse_good_filt = mean((good_eye(:) - filtered_eye(:)).^2);
        mse_difference = mse_good_bad - mse_good_filt;
        %fprintf('Window Size = %d, MSE Difference = %.4f\n', window_size, mse_difference);
        
        % if mse_difference < 0
        %     break;
        % end
        
        % Update Best Window if MSE Diff is Higher
        if mse_difference > best_mse_diff
            best_mse_diff = mse_difference;
            best_window_size = window_size;
        end
    end

    fprintf('Best Window Size = %d\n', best_window_size);
    
    % Apply the Wiener filter once more with the best window size
    %[ch1_filtered, ch2_filtered, ch3_filtered] = WienerFilter(badCh1, badCh2, badCh3, outlier_columns, best_window_size);
    %[ch1_filtered, ch2_filtered, ch3_filtered] = WienerFilter2(badCh1, badCh2, badCh3, outlier_columns, best_window_size);
    [ch1_filtered, ch2_filtered, ch3_filtered] = WienerFilter2Adaptive(badCh1, badCh2, badCh3, outlier_columns, best_window_size);
    

    %Get MSE
    filtered_eye = kSpaceToImage(ch1_filtered, ch2_filtered, ch3_filtered);
    mse_good_filt = mean((good_eye(:) - filtered_eye(:)).^2);
    fprintf('MSE between Reference and Filtered Image: %.4f\n', mse_good_filt);

    mse_diff = mse_good_bad - mse_good_filt;
    fprintf('MSE Diff: %.4f\n', mse_diff);

end

%% LOAD IMAGE SLICES AND PERFORM IMAGE FILTERING

    %Choose which slice to load
    slice = 1;
    [badCh1, badCh2, badCh3, goodCh1, goodCh2, goodCh3] = loadSlice(slice);
    fprintf('\nSlice No: %d\n', slice)
    
    %Compute Outlier Columns
    outlier_columns = findOutlierColumn(badCh1, badCh2, badCh3, 1, 5);    
    %outlier_columns = [38, 54,70, 86, 102, 118]; %HARD-CODE VALUES FOR TESTING
    
    disp(outlier_columns);
    
    % Generate Good and Bad Eye Image
    good_eye = kSpaceToImage(goodCh1, goodCh2, goodCh3);
    bad_eye = kSpaceToImage(badCh1, badCh2, badCh3);
    
    % Wiener Filter
    [ch1_filtered, ch2_filtered, ch3_filtered, filtered_eye] = BestWienerFilter(badCh1, badCh2, badCh3, good_eye, outlier_columns);


% OUTPUTS
%close all;
%set range of colorbar
caxis_range = [0, max(log(abs(badCh1(:))))];

% K-Space Figures
figure;

subplot(3, 1, 1);
imagesc(log(abs(badCh1)));
colorbar;
clim(caxis_range);
title('Channel 1 Bad Data K-Space');
xlabel('Horizontal frequency bins');
ylabel('Vertical frequency bins');

subplot(3, 1, 2);
imagesc(log(abs(ch1_filtered)));
colorbar;
clim(caxis_range);
title('Channel 1 Filtered Data K-Space');
xlabel('Horizontal frequency bins');
ylabel('Vertical frequency bins');

subplot(3, 1, 3);
imagesc(log(abs(goodCh1)));
colorbar;
clim(caxis_range);
title('Channel 1 Good Data K-Space');
xlabel('Horizontal frequency bins');
ylabel('Vertical frequency bins');

% Eye Figures
figure;

subplot(1, 3, 1);
imagesc(bad_eye(:,:,1));
title('Bad Data Eye Image');
axis image;
colormap gray;
axis off;

subplot(1, 3, 2);
imagesc(filtered_eye(:,:,1));
title('Filtered Data Eye Image');
axis image;
colormap gray;
axis off;

subplot(1, 3, 3);
imagesc(good_eye(:,:,1));
title('Good Data Eye Image');
axis image;
colormap gray;
axis off;