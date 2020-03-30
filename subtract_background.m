% Use the the fact that background (mostly small blood vessels) overlaps in different
% channels to subtract it off. 

% Specify one tif file for each channel (filenames should be without
% extension)
folder = {'F:\HCR_9.2_S6'};
filenames = {'HCR_9.2_S6_DAPI';
                'HCR_9.2_S6_Npy2r'; ...
                'HCR_9.2_S6_Tac2'; ...
                'HCR_9.2_S6_Vip'};
            
% Choose which tile to process       
frames = [1];

for frame = frames
    
    % Load data
    raw_channels = load_tifs(folder, filenames, frame);
    n_channel = size(raw_channels, 3);
    

    % Get user to specify background region, mean subtract each channel and normalize to background value

    bg_channel = 2;
    [norm_channels, bg_val, bg_mask_full] = get_bg(raw_channels, n_channel, bg_channel);

    % Choose which channels to use for background calculation
    corr_channels = [2, 3, 4];


    % Compute background as correlation value for each pixel
    plot = 1;
    [v, corr] = get_corr(corr_channels, norm_channels, 'variance', plot);

    % Subtract background and rescale to original image
    sub_channels = get_sub(corr, bg_mask_full, raw_channels, bg_val, folder, frame, n_channel);
    
    % Plot raw and subtracted channels
    figure
    axes = [];
    for i = 1:n_channel
        axes = [axes subplot(n_channel, 2, 2*(i - 1) + 1)];
        imagesc(raw_channels(:, :, i))
        colorbar()
        title(sprintf('Channel %d raw', i))
        
        axes = [axes subplot(n_channel, 2, 2*i)];
        imagesc(sub_channels(:, :, i))
        colorbar()
        title(sprintf('Channel %d background subtracted', i))
    end
    
    linkaxes(axes)

end



