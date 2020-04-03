
from os.path import sep
import pickle as pkl

def get_metadata(data_path, metadata_file, um_per_px, raw_image_path, seg_image_path, plots_path,
                                    filt_image_path, plane_nos, base_filename, genes, channel_names, n_planes, h, w):

    try:
        with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
            metadata = pkl.load(f)
    except:
        metadata = {}

    cell_data_file = 'cell_data.pkl'
    cell_pixels_file = 'cell_pixels.pkl'
    all_cell_pixels_file = 'all_cell_pixels.pkl'
    bg_pixels_file = 'bg_pixels.pkl'
    all_bg_pixels_file = 'all_bg_pixels.pkl'
    lipo_rois_file = 'lipo_rois_cells.pkl'
    lipo_pixels_cells_file = 'lipo_pixels_cells.pkl'

    cell_pixels_filt_file = 'cell_filt_values.pkl'
    bg_pixels_filt_file = 'bg_filt_values.pkl'
    lipo_pixels_cell_filt_file = 'anno_lipo_filt_values.pkl'
    m_dist_cells_file = 'm_dist_cells.pkl'
    m_dist_bg_file = 'm_dist_bg.pkl'
    m_dist_lipo_file = 'm_dist_lipo.pkl'

    cell_pixels_no_lipo_file = 'cell_pixels_no_lipo.pkl'
    bg_pixels_no_lipo_file = 'bg_pixels_no_lipo.pkl'

    depths_file = 'depths.pkl'

    metadata['raw_image_path'] = raw_image_path
    metadata['filt_image_path'] = filt_image_path
    metadata['plane_nos'] = plane_nos
    metadata['base_filename'] = base_filename
    metadata['genes'] = genes
    metadata['channel_names'] = channel_names
    metadata['n_planes'] = n_planes
    metadata['h'] = h
    metadata['w'] = w
    metadata['um_per_px'] = um_per_px

    metadata['cell_data_file'] = cell_data_file
    metadata['cell_pixels_file'] = cell_pixels_file
    metadata['bg_pixels_file'] = bg_pixels_file
    metadata['all_bg_pixels_file'] = all_bg_pixels_file
    metadata['all_cell_pixels_file'] = all_cell_pixels_file
    metadata['lipo_rois_file'] = lipo_rois_file
    metadata['lipo_pixels_cells_file'] = lipo_pixels_cells_file

    metadata['cell_pixels_filt_file'] = cell_pixels_filt_file
    metadata['bg_pixels_filt_file'] = bg_pixels_filt_file
    metadata['lipo_pixels_cells_filt_file'] = lipo_pixels_cell_filt_file
    metadata['m_dist_cells_file'] = m_dist_cells_file
    metadata['m_dist_bg_file'] = m_dist_bg_file
    metadata['m_dist_lipo_file'] = m_dist_lipo_file

    metadata['cell_pixels_no_lipo_file'] = cell_pixels_file
    metadata['bg_pixels_no_lipo_file'] = bg_pixels_file

    metadata['signal_raw_file'] = 'signal_raw.pkl'
    metadata['signal_filt_file'] = 'signal_filt.pkl'
    metadata['bg_raw_file'] = 'bg_raw.pkl'
    metadata['bg_filt_file'] = 'bg_filt.pkl'
    metadata['signal_sigma_small'] = 1
    metadata['signal_sigma_large'] = 2
    metadata['signal_thresh_scale'] = 12

    metadata['seg_image_path'] = seg_image_path
    metadata['plots_path'] = plots_path

    metadata['lipo_sigma_small'] = 5
    metadata['lipo_sigma_large'] = 10
    metadata['lipo_thresh_scale'] = 0.0

    metadata['depths_file'] = depths_file
    metadata['centers_file'] = 'centers.pkl'

    metadata['bs_raw_file'] = 'bs_raw.pkl'
    metadata['bs_filt_file'] = 'bs_filt.pkl'
    metadata['pv_raw_file'] = 'pv_raw.pkl'
    metadata['pv_filt_file'] = 'pv_filt.pkl'


    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'wb') as f:
        pkl.dump(metadata, f)

    return metadata
