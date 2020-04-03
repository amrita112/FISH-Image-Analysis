import matplotlib.pyplot as plt
import numpy as np

def view_all_cells(n_cells, n_rows, n_cols, signal, im_array, cell_data, size_um, um_per_px, h, w,
                  save_loc, title, max_intensity = None, spine_color_pos = 'r', pos_cells = [], spine_width = 2,
                  show_masks = False):

    # Load metadata
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)
    raw_image_path = metadata['raw_image_path']
    filt_image_path = metadata['filt_image_path']
    um_per_px = metadata['um_per_px']
    n_planes = metadata['n_planes']
    h = metadata['h']
    w = metadata['w']
    base_filename = metadata['base_filename']
    channel_names = metadata['channel_names']
    if genes == None:
        genes = metadata['genes']

    # Load cell pixels
    all_cell_pixels_file = metadata['all_cell_pixels_file']
    with open('{0}{1}{2}'.format(data_path, sep, all_cell_pixels_file), 'rb') as f:
        all_cell_pixels = pkl.load(f)

    cell_pixels_no_lipo_file = metadata['cell_pixels_no_lipo_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_no_lipo_file), 'rb') as f:
        cell_pixels_no_lipo = pkl.load(f)

    # Load background pixels
    bg_pixels_no_lipo_file = metadata['bg_pixels_no_lipo_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_pixels_no_lipo_file), 'rb') as f:
        bg_pixels_no_lipo = pkl.load(f)

    # Load signal
    signal_raw_file = metadata['signal_raw_file']
    with open('{0}{1}{2}'.format(data_path, sep, signal_raw_file), 'rb') as f:
        signal_raw = pkl.load(f)

    signal_filt_file = metadata['signal_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, signal_filt_file), 'rb') as f:
        signal_filt = pkl.load(f)

    # Load background
    bg_raw_file = metadata['bg_raw_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_raw_file), 'rb') as f:
        bg_raw = pkl.load(f)

    bg_filt_file = metadata['bg_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_filt_file), 'rb') as f:
        bg_filt = pkl.load(f)

    cells = np.array(list(cell_data.keys()))

    for gene in genes:

        sig = signal_raw[gene]
        sig = np.array([np.mean([sig[cell][plane] for plane in sig[cell].keys()]) for cell in sig.keys()])
        bg = bg_raw[gene]
        bg = np.array([np.mean([bg[cell][plane] for plane in bg[cell].keys()]) for cell in bg.keys()])
        cell_order = cells[np.flip(np.argsort(sig - bg)).astype(int)]

    size_px = int(size_um/um_per_px/2)
    mask_dil =  np.zeros([h, w]).astype(bool)

    im_array = np.zeros(n_planes, h, w])

    for plane in range(n_planes):
        print('Plane {0}, {1} seconds'.format(plane + 1, np.round(time.time() - t0)))
        img = Image.open('{0}{1}{2}{3}{4}.tif'.format(raw_image_path, sep, base_filename, str(plane + 1).zfill(2),
                                                        channel_names[gene]))
        im_array[plane, :, :] = np.array(img)

    if max_intensity == None:
        max_intensity = np.max(sig)
    min_intensity = np.min(im_array)

    t0 = time.time()
    fig, ax = plt.subplots(nrows = n_rows, ncols = n_cols, figsize = (30, 20))

    for cell in cells:

        if np.mod(cell, 10) == 0:
            print('Cell {0}: {1} seconds'.format(cell, np.round(time.time() - t0)))

        row = int(cell/n_cols)
        col = np.mod(cell, n_cols)

        cell_id = cell_order[cell]

        planes = cell_data[cell_id]['z_planes']
        med_plane = int(np.median(planes))

        mask = cell_data[cell_id]['masks'][med_plane]
        [x_center, y_center] = np.mean(mask, 0).astype(int)
        [x1, y1] = [np.max([0, x_center - size_px]), np.max([0, y_center - size_px])]
        [x2, y2] = [np.min([h, x_center + size_px]), np.min([w, y_center + size_px])]
        mask_dil[x1:x2, y1:y2] = np.ones([x2 - x1, y2 - y1])

        small_im = np.reshape(im_array[med_plane, mask_dil], [np.min([x2 - x1, size_px*2]), np.min([y2 - y1, size_px*2])])
        ax[row, col].imshow(small_im, vmin = min_intensity, vmax = max_intensity)

        mask_dil[x1:x2, y1:y2] = np.zeros([x2 - x1, y2 - y1])

        if show_masks:

            cell_px = cell_pixels_no_lipo[cell_id][med_plane]
            cell_mask = np.zeros(small_im.shape)
            x = cell_px[0, :] - x1
            y = cell_px[1, :] - y1
#             for i in range(len(x)):
#                 cell_mask[x[i], y[i]] = 1
            cell_mask[x, y] = np.ones(len(x))
            ax[row, col].imshow(cell_mask, vmin = 0, vmax = 1, cmap = 'gray', alpha = 0.2)


#         ax[row, col].axis('off')
        if col == 0:
            ax[row, col].set_ylabel('{0}'.format(np.round(signal[cell], 2)))
#             ax[row, col].axis('on')

        if len(pos_cells) > 0:
            if pos_cells[cell] == 1:
                ax[row, col].axis('on')
                ax[row, col].spines['bottom'].set_color(spine_color_pos)
                ax[row, col].spines['top'].set_color(spine_color_pos)
                ax[row, col].spines['left'].set_color(spine_color_pos)
                ax[row, col].spines['right'].set_color(spine_color_pos)
                ax[row, col].spines['bottom'].set_linewidth(spine_width)
                ax[row, col].spines['top'].set_linewidth(spine_width)
                ax[row, col].spines['left'].set_linewidth(spine_width)
                ax[row, col].spines['right'].set_linewidth(spine_width)

        ax[row, col].tick_params(axis = 'both', which = 'both', bottom = False, top = False, left = False,
                                 right = False, labeltop = False, labelright = False, labelleft = False,
                                labelbottom = False)


#         plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom=False,      # ticks along the bottom edge are off
#     top=False,         # ticks along the top edge are off
#     labelbottom=False) # labels along the bottom edge are off


    plt.savefig('{0}{1}{2}.png'.format(save_loc, sep, title))
