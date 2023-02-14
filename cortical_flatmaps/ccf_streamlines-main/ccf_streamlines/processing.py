import numpy as np

def remove_duplicate_voxels_from_paths(paths):
    """ Remove duplicate consecutive voxels from the paths

    Duplicate consecutive voxels are identified and removed - the list of
    voxels is then shifted to fill in the gaps left by the duplicates, so the
    format is consistent with the original.

    Parameters
    ----------
    paths : array
        Set of paths containing some duplicated voxels

    Returns
    -------
    new_paths : array
        Set of paths with duplicate voxels removed
    """

    # Remove duplicate consecutive voxels from the paths:
    # First identify places where voxels change (i.e. no duplicates) as
    # values to keep
    paths_diff = np.diff(paths, axis=1)

    # Determine how many remaining entries per line and how many zeros to
    # add at the end
    nonzero_per_row = np.count_nonzero(paths_diff, axis=1)
    n_paths, max_len = paths.shape
    zeros_at_end = max_len - nonzero_per_row

    # Insert the zeros at the ends of each line, then put back together
    # into 2D array
    insert_locs = np.cumsum(nonzero_per_row)
    paths = np.insert(
        paths[np.nonzero(paths_diff)],
        np.repeat(insert_locs, zeros_at_end),
        0).reshape(n_paths, -1)

    return paths