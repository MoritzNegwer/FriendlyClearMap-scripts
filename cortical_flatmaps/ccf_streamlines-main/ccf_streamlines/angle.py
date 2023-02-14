import h5py
import numpy as np
from ccf_streamlines.coordinates import coordinates_to_voxels
from ccf_streamlines.projection import _matching_voxel_indices
from scipy.spatial.distance import euclidean


def vector_to_3d_affine_matrix(vec):
    M = np.array([[vec[0], vec[1], vec[2], vec[9]],
              [vec[3], vec[4], vec[5], vec[10]],
              [vec[6], vec[7], vec[8], vec[11]],
             ])
    return M


def find_closest_streamline(
    coord,
    closest_surface_voxel_reference,
    surface_paths,
    resolution=(10, 10, 10),
    volume_shape=(1320, 800, 1140),
    ):
    """ Find the nearest streamline of a CCF coordinate.

    If voxel is not in the reference file, returns None.

    Parameters
    ----------
    coord : (3, ) array
        Coordinates of point (microns)
    closest_surface_voxel_reference : str or array
        Either a file path to a HDF5 file containing information about the closest
        streamlines for voxels within the isocortex, or the already loaded
        array from that file.
    surface_paths : str or h5py file object
        Either a file path to an HDF5 file containing information about the paths between
        the top and bottom of cortex or an open h5py file object referring to that information.
    resolution : tuple, default (10, 10, 10)
        3-tuple of voxel dimensions in microns
    volume_shape : tuple, default (1320, 800, 1140)
        3-tuple of volume shape in voxels

    Returns
    -------
    streamline_coords : (N, 3) array
        3-D coordinates of streamline in microns
    """

    coord = np.array(coord)
    if len(coord.shape) == 1:
        coord = coord.reshape(1, -1)

    if isinstance(closest_surface_voxel_reference, str):
        with h5py.File(closest_surface_voxel_reference, "r") as f:
            closest_dset = f["closest surface voxel"]
            closest_surface_voxels = closest_dset[:]
    else:
        closest_surface_voxels = closest_surface_voxel_reference

    voxel = np.squeeze(coordinates_to_voxels(coord))

    # Reference file data only present on left side, so flip to left side
    # if voxel is on the right
    z_size = volume_shape[2]
    z_midline = z_size / 2
    flip_hemisphere = False
    if voxel[2] > z_midline:
        voxel[2] = z_size - voxel[2]
        flip_hemisphere = True

    voxel_ind = np.ravel_multi_index(
            tuple(voxel),
            volume_shape
        )
    matching_surface_voxel_ind = _matching_voxel_indices(
        np.array([voxel_ind]),
        closest_surface_voxels)[0]

    # Pull path from surface paths file
    if isinstance(surface_paths, str):
        with h5py.File(surface_paths, "r") as f:
            path_dset = f['paths']
            volume_lookup_dset = f['volume lookup flat']
            path_ind = volume_lookup_dset[matching_surface_voxel_ind]
            path = path_dset[path_ind, :]
    else:
        path_dset = surface_paths['paths']
        volume_lookup_dset = surface_paths['volume lookup flat']
        path_ind = volume_lookup_dset[matching_surface_voxel_ind]
        path = path_dset[path_ind, :]

    path = path[path > 0]

    # Convert path to coordinates in microns
    streamline_coords = np.unravel_index(
        path, volume_shape)
    streamline_coords = np.array(streamline_coords).T
    if flip_hemisphere:
        # Put streamline on same hemisphere as original voxel
        streamline_coords[:, 2] = z_size - streamline_coords[:, 2]

    # Scale to microns
    streamline_coords = streamline_coords * np.array(resolution)
    return streamline_coords


def determine_angle_between_streamline_and_plane(
    streamline_coords,
    plane_transform):
    """
    Find angle between the surface of a given plane and a streamline.

    Parameters
    ----------
    streamline_coords : (N, 3) array
        Coordinates in microns of streamline
    plane_transform : (3, 4) array
        3-D affine transform matrix for putting a plane into CCF space.

    Returns
    -------
    angle : float
        Angle between streamline and plane (in degrees)
    """

    # Find normal vector for plane transformed to CCF space
    a2 = np.dot(plane_transform, np.array([0, 0, 0, 1]))
    b2 = np.dot(plane_transform, np.array([1, 0, 0, 1]))
    c2 = np.dot(plane_transform, np.array([0, 1, 0, 1]))

    norm_vec = np.cross(b2 - a2, c2 - a2)
    norm_unit = norm_vec / euclidean(norm_vec, [0, 0, 0])

    streamline_wm = streamline_coords[-1, :]
    streamline_pia = streamline_coords[0, :]
    streamline_unit = (streamline_pia - streamline_wm) / euclidean(streamline_pia, streamline_wm)

    angle_with_norm = np.arctan2(
        np.linalg.norm(np.cross(norm_unit, streamline_unit)),
        np.dot(norm_unit, streamline_unit))

    # Take complement and convert from radians to degrees
    return 90. - (angle_with_norm * 180. / np.pi)
