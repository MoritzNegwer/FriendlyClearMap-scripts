import numpy as np
import pandas as pd
import nrrd
import h5py
import logging
from tqdm import tqdm


def measure_streamline_layer_thicknesses(layer_volume, paths, resolution):
    """ Measure the start, end, and thickness of layers

    Parameters
    ----------
    layer_volume : array
        3D volume with annotated layers
    paths : array
        Streamline paths
    resolution : tuple
        3-tuple of voxel size in x, y, z

    Returns
    -------
    thicknesses : dict
        Dictionary keyed on layers with start, end, and thickness of layers
    """

    # Structure set IDs from mouse ontology
    LAYER_LABELS = {
        667481440: 'Isocortex layer 1',
        667481441: 'Isocortex layer 2/3',
        667481445: 'Isocortex layer 4',
        667481446: 'Isocortex layer 5',
        667481449: 'Isocortex layer 6a',
        667481450: 'Isocortex layer 6b',
    }

    # Remove duplicate consecutive voxels from paths
    fixed_paths = np.zeros_like(paths)
    paths_diff = np.diff(paths, axis=1)
    for i in range(paths.shape[0]):
        unique_inds = np.flatnonzero(paths_diff[i, :])
        fixed_paths[i, :len(unique_inds)] = paths[i, :][unique_inds]
    paths = fixed_paths

    max_nonzero_path_inds = np.count_nonzero(paths, axis=1)

    # get voxel coordinates for all path voxels
    path_x, path_y, path_z = np.unravel_index(paths, layer_volume.shape)
    path_voxels = np.stack([
        path_x * resolution[0],
        path_y * resolution[1],
        path_z * resolution[2]
    ])

    # Find distances between consecutive voxels
    deltas = np.diff(path_voxels, axis=2)
    distances = np.sqrt((deltas ** 2).sum(axis=0))

    # add thickness of last voxel to end of paths
    distances[np.arange(distances.shape[0]), max_nonzero_path_inds - 1] = np.mean(resolution)

    # Calculate cumulative depths
    cumul_distances = np.cumsum(distances, axis=1)

    layer_annot = layer_volume.flat[paths]

    # Sort the layer annotations so that layers cannot be intercalated
    # (Obviously this is a real hack, but don't think there's a robust
    # solution to the issue)
    # Set zeros to large value so they end up last
    layer_annot[layer_annot == 0] = 999999999
    layer_annot = np.sort(layer_annot, axis=1)

    thicknesses = {}
    end_depths = {}
    start_depths = {}
    for k, v in LAYER_LABELS.items():
        layer_mask = layer_annot == k

        # Last row dropped from np.diff
        layer_mask = layer_mask[:, :-1]

        thicknesses[v] = np.sum(np.where(layer_mask, distances, 0), axis=1)
        end_depths[v] = np.max(np.where(layer_mask, cumul_distances, 0), axis=1)
        start_depths[v] = end_depths[v] - thicknesses[v]

    output = {}
    for k in thicknesses:
        output[k] = np.vstack([start_depths[k], end_depths[k], thicknesses[k]]).T

        # set very small values to exactly zero
        output[k][np.isclose(output[k], 0)] = 0
    return output
