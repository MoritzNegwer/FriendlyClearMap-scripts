import numpy as np
import pandas as pd
from ccf_streamlines.coordinates import coordinates_to_voxels


def load_swc_as_dataframe(swc_file):
    """ Load a morphology SWC file into a pandas DataFrame.

    The dataframe contains the columns:
        - ID : node identifier
        - type : node type, which could be
            1 = soma
            2 = axon
            3 = basal dendrite (or generic dendrite)
            4 = apical dendrite
        - x : x-coordinate
        - y : y-coordinate
        - z : z-coordinate
        - r : radius
        - parent_id : identifier of parent node

    Parameters
    ----------
    swc_file : str
        File path of SWC file

    Returns
    -------
    df : DataFrame
        Dataframe with morphology information
    """
    return pd.read_table(
        swc_file,
        sep=" ",
        comment="#",
        names=["id", "type", "x", "y", "z", "r", "parent_id"],
    )


def transform_swc_to_volume(
    swc_file,
    volume_shape=(1320, 800, 1140),
    resolution=(10, 10, 10),
    compartments=None
    ):
    """ Create a volume with node counts from a morphology SWC file.

    Parameters
    ----------
    swc_file : str
        File path of SWC file
    volume_shape : 3-tuple of ints, default (1320, 800, 1140)
        Shape of target volume in voxels
    resolution : 3-tuple of ints, default (10, 10, 10)
        Size of voxels in microns
    compartments : list-like, optional
        List of compartment types to include
        1 = soma
        2 = axon
        3 = basal dendrite (or generic dendrite)
        4 = apical dendrite
        If None (default), all compartments are included.

    Returns
    -------
    volume : array
        Volume with shape `volume_shape` with counts of nodes in each voxel
    """

    # Load the coordinates from the SWC file
    swc_df = load_swc_as_dataframe(swc_file)

    # Limit to requested compartment types
    if compartments is not None:
        swc_df = swc_df.loc[swc_df['type'].isin(compartments), :]

    # Determine voxels from coordinates
    coords = swc_df.loc[:, ['x', 'y', 'z']].values

    return transform_coordinates_to_volume(coords, volume_shape, resolution)


def transform_coordinates_to_volume(
    coords,
    volume_shape=(1320, 800, 1140),
    resolution=(10, 10, 10),
    ):
    """ Create a volume with counts of points.

    Parameters
    ----------
    coords : (N, 3) array
        Coordinates of points
    volume_shape : 3-tuple of ints, default (1320, 800, 1140)
        Shape of target volume in voxels
    resolution : 3-tuple of ints, default (10, 10, 10)
        Size of voxels in microns

    Returns
    -------
    volume : array
        Volume with shape `volume_shape` with counts of points in each voxel
    """
    voxels = coordinates_to_voxels(coords, resolution=resolution)

    # Count the nodes in each voxel
    populated_voxels, counts = np.unique(voxels, axis=0, return_counts=True)

    # Place counts into volume
    volume = np.zeros(volume_shape, dtype=np.uint32)
    volume[populated_voxels[:, 0], populated_voxels[:, 1], populated_voxels[:, 2]] = counts

    return volume


def find_topological_point_coordinates(swc_df):
    """ Get the coordinates of the branch and termination nodes

    Parameters
    ----------
    swc_df : dataframe
        Dataframe with morphology information

    Returns
    -------
    coords : (N, 3) array
        3-D coordinates of branch and termination nodes
    """

    child_counts = swc_df["parent_id"].value_counts()
    branch_ids = child_counts[child_counts > 1].index.intersection(swc_df["id"]).tolist()
    term_ids = [i for i in swc_df["id"] if i not in child_counts.index]

    return swc_df.set_index("id").loc[branch_ids + term_ids, ["x", "y", "z"]].values


