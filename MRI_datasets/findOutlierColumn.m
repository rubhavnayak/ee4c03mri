%Function to detect outlier columns by analyzing MAD of three channels
%Arguments: Three channels of K-space images, peak height multiplier, minimum peak distance
function outlier_columns = findOutlierColumn(ch1, ch2, ch3, n_std, minPeakDistance)
   
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
        [~, locs] = findpeaks(detrended_median_diff, 'MinPeakHeight', n_std*std(detrended_median_diff), 'MinPeakDistance', minPeakDistance);
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