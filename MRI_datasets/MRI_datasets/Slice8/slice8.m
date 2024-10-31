%%CLEAR OUTPUTS AND WORKSPACE
clc; clear all; close all;

%% Load K-space data for each channel
badCh1 = load('BadData/slice8_channel1.mat', 'slice8_channel1_badData').slice8_channel1_badData;
badCh2 = load('BadData/slice8_channel2.mat', 'slice8_channel2_badData').slice8_channel2_badData;
badCh3 = load('BadData/slice8_channel3.mat', 'slice8_channel3_badData').slice8_channel3_badData;

goodCh1 = load('GoodData/slice8_channel1.mat', 'slice8_channel1_goodData').slice8_channel1_goodData;
goodCh2 = load('GoodData/slice8_channel2.mat', 'slice8_channel2_goodData').slice8_channel2_goodData;
goodCh3 = load('GoodData/slice8_channel3.mat', 'slice8_channel3_goodData').slice8_channel3_goodData;

%% INITIALIZE FUNCTIONS

%Function to fuse 3 channels of K-space images to obtain reconstructed eye image.
%Arguments: Three channels of K-space images.
function eye_visualize = kSpaceToImage(ch1, ch2, ch3)
    
    % IFFT of k-space data
    Data_img(:,:,1) = ifftshift(ifft2(ch1),1); %channel 1
    Data_img(:,:,2) = ifftshift(ifft2(ch2),1); %channel 2
    Data_img(:,:,3) = ifftshift(ifft2(ch3),1); %channel 3
    
    % clear compensation, preparation, based on fourier transformed blinked k-space data (Data_raw)
    clear_comp = linspace(10,0.1,size(Data_img,2)).^2; 
    clear_matrix = repmat(clear_comp,[size(Data_img,1) 1]);
    
    % combine 3 channels sum of squares and add clear compensation
    eye_raw  = sqrt( abs(squeeze(Data_img(:,:,1))).^2 + ... % (.^ signifies element-wise power operation)
               abs(squeeze(Data_img(:,:,2))).^2 + ...
               abs(squeeze(Data_img(:,:,3))).^2).* clear_matrix;  
        
    % crop images because we are only interested in eye. Make it square 
    % 128 x 128
    crop_x = [128 + 60 : 348 - 33]; % crop coordinates
    eye_raw = eye_raw(crop_x, :);
    
    % Visualize the images. 
    eye_visualize = reshape(squeeze(eye_raw(:,:)),[128 128]); 
    
    %Histogram compensation for improved visualization
    std_within = 0.995; % set maximum intensity to contain 99.5 % of intensity values per image
    [aa, val] = hist(eye_visualize(:),linspace(0,max(...
                                        eye_visualize(:)),1000));
    thresh = val(find(cumsum(aa)/sum(aa) > std_within,1,'first'));
        
    % set threshold value to 65536
    eye_visualize = uint16(eye_visualize * 65536 / thresh); 

end

%Function to detect outlier columns by analyzing MAD of three channels
%Arguments: Three channels of K-space images.
function outlier_columns = findOutlierColumn(ch1, ch2, ch3)
   
    outlier_columns_ch1 = [];
    outlier_columns_ch2 = [];
    outlier_columns_ch3 = [];
    
    channels = {ch1, ch2, ch3};
    for ch = 1:3
        channel = channels{ch};
        
        % Compute the median difference of adjacent columns
        diff_combined = abs(diff(channel, 1, 2));
        median_diff = median(diff_combined, 1);

        % Define columns to analyze, excluding central columns
        column_indices = 2:size(channel, 2);
        central_columns = 62:68;
        columns_to_analyze = setdiff(column_indices, central_columns);
        median_diff_to_analyze = median_diff(columns_to_analyze - 1);

        % Detrend data to remove bell-curve trend
        p = polyfit(columns_to_analyze, median_diff_to_analyze, 2);
        trend = polyval(p, columns_to_analyze);
        detrended_median_diff = median_diff_to_analyze - trend;

        % Detect peaks using findpeaks
        [~, locs] = findpeaks(detrended_median_diff, 'MinPeakHeight', std(detrended_median_diff), 'MinPeakDistance', 5);
        outlier_columns = columns_to_analyze(locs);

        % Store corrupted columns for each channel
        if ch == 1
            outlier_columns_ch1 = outlier_columns;
        elseif ch == 2
            outlier_columns_ch2 = outlier_columns;
        else
            outlier_columns_ch3 = outlier_columns;
        end
    end

    % Consistent corrupted columns across channels
    %outlier_columns = intersect(intersect(outlier_columns_ch1, outlier_columns_ch2), outlier_columns_ch3) - 1;
    outlier_columns = outlier_columns - 1;
end

%Median Filter
%Arguments: Three channels of K-space images, vector of outlier columns
function [ch1_filtered, ch2_filtered, ch3_filtered] = medianFilter(ch1, ch2, ch3, outlier_columns)
    % Corrects corrupted columns by replacing them with the median of neighboring columns.
    
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
end

function [ch1_filtered, ch2_filtered, ch3_filtered] = WienerFilter(ch1, ch2, ch3, outlier_columns, window_size)
    % Wiener filter to correct corrupted columns based on local statistics.
    % Inputs:
    % - ch1, ch2, ch3: K-space channels
    % - outlier_columns: vector of corrupted column indices
    % - window_size: size of the neighborhood for local statistics estimation
    % Outputs:
    % - ch1_filtered, ch2_filtered, ch3_filtered: filtered channels

    % Initialize filtered channels
    filtered_channels = {ch1, ch2, ch3};

    for idx = 1:3
        channel = filtered_channels{idx};
        
        % Estimate global noise variance from the entire channel
        noise_variance = var(channel(:));

        for col_idx = outlier_columns
            % Define the neighborhood window for the corrupted column
            window_start = max(col_idx - window_size, 1);
            window_end = min(col_idx + window_size, size(channel, 2));
            window_data = channel(:, window_start:window_end);
            
            % Calculate local mean and variance within the window
            local_mean = mean(window_data, 2);
            local_variance = var(window_data, 0, 2);
            
            % Apply Wiener filter formula
            estimated_column = local_mean + ...
                (max(local_variance - noise_variance, 0) ./ max(local_variance, noise_variance)) .* ...
                (window_data(:, ceil(end/2)) - local_mean);
            
            % Replace corrupted column with Wiener-filtered data
            channel(:, col_idx) = estimated_column;
        end
        
        filtered_channels{idx} = channel;
    end

    % Output filtered channels
    ch1_filtered = filtered_channels{1};
    ch2_filtered = filtered_channels{2};
    ch3_filtered = filtered_channels{3};
end

%% FILTER THE IMAGES AND GET THE OUTPUTS

%Compute Outlier Columns
outlier_columns = findOutlierColumn(badCh1, badCh2, badCh3);

%Hard-coded Outlier Columns (gives a slightly better result)
%outlier_columns = [38, 54,70, 86, 102, 118]; %HARD-CODE VALUES

%Choice of Filters
%[ch1_filtered, ch2_filtered, ch3_filtered] = medianFilter(badCh1, badCh2, badCh3, outlier_columns);
[ch1_filtered, ch2_filtered, ch3_filtered] = WienerFilter(badCh1, badCh2, badCh3, outlier_columns, 59);

% OUTPUTS
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

bad_eye = kSpaceToImage(badCh1, badCh2, badCh3);
subplot(1, 3, 1);
imagesc(bad_eye(:,:,1));
title('Bad Data Eye Image');
axis image;
colormap gray;
axis off;

filtered_eye = kSpaceToImage(ch1_filtered, ch2_filtered, ch3_filtered);
subplot(1, 3, 2);
imagesc(filtered_eye(:,:,1));
title('Filtered Data Eye Image');
axis image;
colormap gray;
axis off;

good_eye = kSpaceToImage(goodCh1, goodCh2, goodCh3);
subplot(1, 3, 3);
imagesc(good_eye(:,:,1));
title('Good Data Eye Image');
axis image;
colormap gray;
axis off;

% MSE CALCULATIONS
mse_good_bad = mean((good_eye(:) - bad_eye(:)).^2);
fprintf('MSE between Reference and Original Image: %.4f\n', mse_good_bad);

mse_good_filt = mean((good_eye(:) - filtered_eye(:)).^2);
fprintf('MSE between Reference and Filtered Image: %.4f\n', mse_good_filt);

mse_diff = mse_good_bad - mse_good_filt;
fprintf('MSE Diff: %.4f\n', mse_diff);