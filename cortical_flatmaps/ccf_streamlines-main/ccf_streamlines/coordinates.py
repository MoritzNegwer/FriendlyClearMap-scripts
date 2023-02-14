import numpy as np


def coordinates_to_voxels(coords, resolution=(10, 10, 10)):
    """ Find the voxel coordinates of spatial coordinates

    Parameters
    ----------
    coords : array
        (n, m) coordinate array. m must match the length of `resolution`
    resolution : tuple, default (10, 10, 10)
        Size of voxels in each dimension

    Returns
    -------
    voxels : array
        Integer voxel coordinates corresponding to `coords`
    """

    if len(resolution) != coords.shape[1]:
        raise ValueError(
            f"second dimension of `coords` must match length of `resolution`; "
            f"{len(resolution)} != {coords.shape[1]}")

    if not np.issubdtype(coords.dtype, np.number):
        raise ValueError(f"coords must have a numeric dtype (dtype is '{coords.dtype}')")

    voxels = np.floor(coords / resolution).astype(int)
    return voxels
