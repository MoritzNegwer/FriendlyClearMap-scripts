import numpy as np
import pytest
import ccf_streamlines.coordinates as coordinates


def test_mismatch():
    test_coords = np.array([
        [0, 0, 0, 0],
        [1, 1, 1, 1],
    ])
    resolution = (10, 10, 10)

    with pytest.raises(ValueError):
        coordinates.coordinates_to_voxels(
            test_coords,
            resolution)


def test_coords_to_voxels():
    test_coords = np.array([
        [0., 0., 0.],
        [5., 5., 5.],
        [10., 10., 10.],
        [15., 25., 0.],
    ])
    expected_voxels = np.array([
        [0, 0, 0],
        [0, 0, 0],
        [1, 1, 1],
        [1, 2, 0],
    ])

    resolution = (10, 10, 10)

    assert np.all(
        coordinates.coordinates_to_voxels(
            test_coords,
            resolution
        ) == expected_voxels
    )

    # Try doubling coordinates and resolution;
    # should get same answer
    double_resolution = (20, 20, 20)
    assert np.all(
        coordinates.coordinates_to_voxels(
            test_coords * 2,
            double_resolution
        ) == expected_voxels
    )
