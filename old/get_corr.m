    
function [v, corr] = get_corr(corr_channels, channels, corr_type, plot)

    if nargin < 4
        plot = 1;
    end
    
    if nargin < 3
        corr_type = 'variance';
    end

    corr = ones(size(channels, 1), size(channels, 2), 1);
    
    if strcmp(corr_type, 'product')
        for i = corr_channels
            corr = corr.*channels(:, :, i);
        end

        corr(corr < 0) = 0;
        v = zeros(size(corr));
    end
    
    if strcmp(corr_type, 'variance')
        channels2 = zeros(size(channels, 3), size(channels, 1), size(channels, 2));
        for i = 1:size(channels, 3)
            channels2(i, :, :) = channels(:, :, i);
        end
        v = var(channels2);
        %v = v(v > thresh);
        v = squeeze(v);
        corr = ones(size(v))./v;
    end

        
    if plot
        figure
        imagesc(corr)
        title('Correlation across channels')
        colorbar()
    end
end