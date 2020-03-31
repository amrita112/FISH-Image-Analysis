
function [channels, bg_val, bg_mask_full] = get_bg(channels, n_channel, bg_channel)


    % Ask user to select a region with autofluorescence and draw a ROI
    % around it
    figure
    imagesc(channels(:, :, bg_channel))
    title('Draw a rectangle around a small area with autofluorescence')
    [im_bg, rect2] = imcrop();
    close 
    
    figure
    imagesc(im_bg)
    title('Draw a polygon over the brightest background region')
    bg_mask = roipoly();
    
    % Extend the mask to the full image
    bg_mask_full = zeros(size(channels, 1), size(channels, 2));
    x1 = round(rect2(1));
    x2 = x1 + round(rect2(3));
    y1 = round(rect2(2));
    y2 = y1 + round(rect2(4));
    bg_mask_full(y1:y2 - 1, x1:x2 - 1) = bg_mask;
    bg_mask_full = logical(bg_mask_full);
    
    % Get background value for each channel
    bg_val = zeros(n_channel, 1);
    mean_val = zeros(n_channel, 1);

    for i = 1:n_channel
        mean_val(i) = mean(mean(channels(:, :, i)));
        channels(:, :, i) = channels(:, :, i) - mean_val(i);
        img = channels(:, :, i);
        bg = img(bg_mask_full);
        bg_val(i) = mean(bg);
        disp(strcat('Background for channel ', int2str(i), ' = ', string(bg_val(i))))
        channels(:, :, i) = channels(:, :, i)/bg_val(i);
    end
end