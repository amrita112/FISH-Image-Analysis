% One tile per channel only


function channels = load_tifs(folder, filenames, frame)
    
    n_channel = numel(filenames);
%    folder = {'F:\452392'};
    % Filenames should be without extension
 %   filenames = {'S3_594';
  %                  'S3_647';
   %                 'S3_jf525'};

    image = Tiff(strjoin([folder '\' string(filenames(1)) '.tif'], ''), 'r');
    im_size = [image.getTag('ImageLength'), image.getTag('ImageWidth')];
    channels = zeros(im_size(1), im_size(2), n_channel);

    for i = 1:n_channel
            disp(i)
            %channels(:, :, i) = read(Tiff(strjoin([folder '\' string(filenames(i)) '.tif'], ''), 'r'));
            channels(:, :, i) = imread(strjoin([folder '\' string(filenames(i)) '.tif'], ''), frame);
    end

    
end