function sub_channels = get_sub(corr, bg_mask, raw_channels, bg_val, folder, frame, n_channel)

    bg_corr = mean(corr(bg_mask));
    sub_channels = raw_channels;

    for i = 1:n_channel
        sub_channels(:, :, i) = raw_channels(:, :, i) - bg_val(i)/bg_corr*corr;
        sub_channels(:, :, i) = sub_channels(:, :, i) - min(min(sub_channels(:, :, i)));
        sub_channels(:, :, i) = sub_channels(:, :, i)/max(max(sub_channels(:, :, i)));
    end
    % sub_channels(sub_channels < 0) = 0;
    for i = 1:n_channel

        imwrite(sub_channels(:, :, i), strjoin([folder '\frame_' string(frame) '_channel_' string(i) '_bg_subtract.tif'], ''))
        disp('Background subtracted images saved as:')
        disp(strjoin([folder '\frame_' string(frame) '_channel_' string(i) '_bg_subtract.tif']))

    end

end
