% Function to fuse 3 channels of K-space images to obtain reconstructed eye image.
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
