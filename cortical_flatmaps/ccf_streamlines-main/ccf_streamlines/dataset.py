import numpy as np


def upscale_ish_volume(
    volume,
    orig_voxel_size=200,
    target_voxel_size=10,
    target_volume_shape=(1320, 800, 1140),
    rotate_axes=True,
    ):
    """ Upscale a lower-resolution ISH volume for projection.

    Parameters
    ----------
    volume : array
        Array of ISH expression data
    orig_voxel_size : int, default 200
        Size in microns of voxels of original data
    target_voxel_size : int, default 10
        Size of microns of voxels in desired volume
    target_volume_shape : tuple, default (1320, 800, 1140)
        Shape of target upscaled volume
    rotate_axes : bool, default True
        Whether to swap the x and z axes of the input volume. Volumes downloaded
        directly from the ISH atlas API have the anterior-posterior axis in the
        z-axis and the left-right axis in the x-axis, while the CCFv3 has those
        swapped. The default assumes that the input volume has not been adjusted
        to match the CCF yet.

    Returns
    -------
    upscaled_volume : array
        Upscaled ISH expression data
    """

    if rotate_axes:
        volume = np.swapaxes(volume, 0, 2)

    voxel_ratio = orig_voxel_size // target_voxel_size
    ind = np.indices(target_volume_shape, sparse=True)
    downscaled_ind = tuple(ind_dim // voxel_ratio for ind_dim in ind)
    upscaled_volume = np.zeros(target_volume_shape)
    upscaled_volume[ind] = volume[downscaled_ind]

    return upscaled_volume