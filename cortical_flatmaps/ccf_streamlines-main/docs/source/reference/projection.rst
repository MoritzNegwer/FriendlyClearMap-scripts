.. _projection:

.. module:: ccf_streamlines.projection

Streamline Projection (:mod:`ccf_streamlines.projection`)
=========================================================

This module contains objects that are configured with different flattened projections, 
hemisphere options, and depth normalization functions and then used to project volumes
or sets of coordinates into the flattened representations.

.. currentmodule:: ccf_streamlines

Projecting 3D volumes
---------------------

.. autosummary::
    :toctree: generated/

    projection.Isocortex2dProjector
    projection.Isocortex3dProjector


Projecting sets of coordinates
------------------------------

.. autosummary::
    :toctree: generated/

    projection.IsocortexCoordinateProjector


Identifying cortical area boundaries
------------------------------------

.. autosummary::
    :toctree: generated/

    projection.BoundaryFinder

