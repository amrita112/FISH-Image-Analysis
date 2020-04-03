
from os.path import sep
import numpy as np
import pickle as pkl
from scipy.spatial.distance import mahalanobis
from utils import find_threshold

def get_m_dist(data_path, metadata_file, genes = None):

    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)

    # Get cell pixel values
    cell_pixels_filt_file = metadata['cell_pixels_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_filt_file), 'rb') as f:
        cell_px = pkl.load(f)

    cell_data_file = metadata['cell_data_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_data_file), 'rb') as f:
        cell_data = pkl.load(f)

    # Get background pixel values
    bg_pixels_filt_file = metadata['bg_pixels_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_pixels_filt_file), 'rb') as f:
        bg_px = pkl.load(f)

    # Get annotated lipofuscin pixel values
    lipo_pixels_cells_filt_file = metadata['lipo_pixels_cells_filt_file']
    with open('{0}{1}{2}'.format(data_path, sep, lipo_pixels_cells_filt_file), 'rb') as f:
        lipo_px = pkl.load(f)

    # Concatenate lipofuscin pixel values
    if genes == None:
        genes = list(lipo_px.keys())
    n_genes = len(genes)
    lipo_px_vals = np.zeros([n_genes, len(lipo_px[genes[0]])])
    for g in range(n_genes):
        lipo_px_vals[g, :] = lipo_px[genes[g]]

    mu = np.mean(lipo_px_vals, axis = 1)
    cov = np.cov(lipo_px_vals)
    VI = np.linalg.inv(cov)

    m_dist_cells = {}
    m_dist_bg = {}
    m_dist_lipo = np.zeros(lipo_px_vals.shape[1])

    cells = list(cell_data.keys())
    for cell in cells:

        m_dist_cells[cell] = {}
        m_dist_bg[cell] = {}
        planes = cell_data[cell]['z_planes']
        for plane in planes:
            n_px = len(cell_px[genes[0]][cell][plane])
            u_cell = np.zeros([n_genes, n_px])

            n_px_bg = len(bg_px[genes[0]][cell][plane])
            u_bg = np.zeros([n_genes, n_px_bg])

            for g in range(n_genes):
                u_cell[g, :] = cell_px[genes[g]][cell][plane]
                u_bg[g, :] = bg_px[genes[g]][cell][plane]

            m_dist_cells[cell][plane] = np.zeros(n_px)
            m_dist_bg[cell][plane] = np.zeros(n_px_bg)

            for p in range(n_px):
                m_dist_cells[cell][plane][p] = mahalanobis(u_cell[:, p], mu, VI)
            for p in range(n_px_bg):
                m_dist_bg[cell][plane][p] = mahalanobis(u_bg[:, p], mu, VI)

    n_px = lipo_px_vals.shape[1]
    for p in range(n_px):
        u = lipo_px_vals[:, p]
        m_dist_lipo[p] = mahalanobis(u, mu, VI)

    m_dist_cells_file = metadata['m_dist_cells_file']
    with open('{0}{1}{2}'.format(data_path, sep, m_dist_cells_file), 'wb') as f:
        pkl.dump(m_dist_cells, f)

    m_dist_bg_file = metadata['m_dist_bg_file']
    with open('{0}{1}{2}'.format(data_path, sep, m_dist_bg_file), 'wb') as f:
        pkl.dump(m_dist_bg, f)

    m_dist_lipo_file = metadata['m_dist_lipo_file']
    with open('{0}{1}{2}'.format(data_path, sep, m_dist_lipo_file), 'wb') as f:
        pkl.dump(m_dist_lipo, f)

    return m_dist_cells, m_dist_bg, m_dist_lipo




def get_lipo(data_path, metadata_file, thresh_scale = 1.5):

    # Load metadata
    with open('{0}{1}{2}'.format(data_path, sep, metadata_file), 'rb') as f:
        metadata = pkl.load(f)

    # Load cell mask data
    cell_data_file = metadata['cell_data_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_data_file), 'rb') as f:
        cell_data = pkl.load(f)

    cell_pixels_file = metadata['cell_pixels_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_file), 'rb') as f:
        cell_pixels = pkl.load(f)

    # Load background pixels
    bg_pixels_file = metadata['bg_pixels_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_pixels_file), 'rb') as f:
        bg_pixels = pkl.load(f)

    # Load mahalanobis distances of cell and background pixels
    m_dist_cells_file = metadata['m_dist_cells_file']
    with open('{0}{1}{2}'.format(data_path, sep, m_dist_cells_file), 'rb') as f:
        m_dist_cells = pkl.load(f)

    m_dist_bg_file = metadata['m_dist_bg_file']
    with open('{0}{1}{2}'.format(data_path, sep, m_dist_bg_file), 'rb') as f:
        m_dist_bg = pkl.load(f)

    # Concatenate cell and bg pixel distances to get threshold
    cells = list(cell_data.keys())
    cell = cells[0]
    planes = cell_data[cell]['z_planes']
    points = m_dist_cells[cell][planes[0]]
    points = np.concatenate([points, m_dist_bg[cell][planes[0]]])
    for plane in planes[1:]:
        points = np.concatenate([points, m_dist_cells[cell][plane]])
        points = np.concatenate([points, m_dist_bg[cell][plane]])

    for cell in cells[1:]:
        planes = cell_data[cell]['z_planes']
        for plane in planes:
            points = np.concatenate([points, m_dist_cells[cell][plane]])
            points = np.concatenate([points, m_dist_bg[cell][plane]])

    thresh = find_threshold.find_threshold(points, thresh_scale)

    # Find new cell and background masks excluding Lipofuscin
    cell_pixels_no_lipo = {}
    bg_pixels_no_lipo = {}
    for cell in cells:
        cell_pixels_no_lipo[cell] = {}
        bg_pixels_no_lipo[cell] = {}
        planes = cell_data[cell]['z_planes']
        for plane in planes:
            cell_pixels_no_lipo[cell][plane] = np.array(cell_pixels[cell][plane])[:, np.where(m_dist_cells[cell][plane] > thresh)[0]]
            bg_pixels_no_lipo[cell][plane] = np.array(bg_pixels[cell][plane])[np.where(m_dist_bg[cell][plane] > thresh)[0], :]

    cell_pixels_no_lipo_file = metadata['cell_pixels_no_lipo_file']
    with open('{0}{1}{2}'.format(data_path, sep, cell_pixels_no_lipo_file), 'wb') as f:
        pkl.dump(cell_pixels_no_lipo, f)

    bg_pixels_no_lipo_file = metadata['bg_pixels_no_lipo_file']
    with open('{0}{1}{2}'.format(data_path, sep, bg_pixels_no_lipo_file), 'wb') as f:
        pkl.dump(bg_pixels_no_lipo, f)

    return cell_pixels_no_lipo, bg_pixels_no_lipo
