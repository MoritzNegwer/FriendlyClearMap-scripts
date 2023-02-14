Data Files
==========

The ``ccf_streamlines`` package uses several data files that contain the different flattened
representations, cortical streamlines, and other information.

To use many aspects of ``ccf_streamlines``, these files must be downloaded to the user's computer,
and their file paths are passed as parameters to various objects and functions.

This document links to those files and describes what they contain.


Streamlines
-----------

.. list-table:: Streamline files
    :widths: 25 75
    :header-rows: 1

    * - File
      - Description
    * - `surface_paths_10_v3.h5 <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/streamlines/surface_paths_10_v3.h5>`_
      - Contains the (linear) voxel locations of each streamline and a lookup table to find a streamline for each voxel. **Warning:** Large file (~0.5 GB)
    * - `closest_surface_voxel_lookup.h5 <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/streamlines/closest_surface_voxel_lookup.h5>`_
      - Lookup table of nearest surface voxel for cortical voxels, which are used to find the best-matching streamline. **Warning:** Large file (~200 MB)


View lookups
------------

.. list-table:: View lookup files
    :widths: 25 50 25
    :header-rows: 1

    * - File
      - Description
      - Image
    * - `back.h5 <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/view_lookup/back.h5>`_
      - Correspondence between view from back of isocortex and CCF volume. Can place hemispheres adjacent to each other with ``view_space_for_other_hemisphere = True``.
      - .. image:: /images/back_template.png
           :width: 400
    * - `bottom.h5 <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/view_lookup/bottom.h5>`_
      - Correspondence between view from bottom of isocortex and CCF volume. Can place hemispheres adjacent to each other with ``view_space_for_other_hemisphere = True``.
      - .. image:: /images/bottom_template.png
           :width: 400
    * - `flatmap_butterfly.h5 <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/view_lookup/flatmap_butterfly.h5>`_
      - Correspondence between flattened map of all of isocortex and CCF volume - medial side has been adjusted so that the hemispheres are abutting at the center. Can place hemispheres adjacent to each other with ``view_space_for_other_hemisphere = 'flatmap_butterfly'``.
      - .. image:: /images/flatmap_butterfly_template.png
           :width: 400
    * - `flatmap_dorsal.h5 <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/view_lookup/flatmap_dorsal.h5>`_
      - Correspondence between flattened map of all of isocortex and CCF volume. Can place hemispheres near each other with ``view_space_for_other_hemisphere = 'flatmap_dorsal'``.
      - .. image:: /images/flatmap_dorsal_template.png
           :width: 400
    * - `front.h5 <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/view_lookup/front.h5>`_
      - Correspondence between view from front of isocortex and CCF volume. Can place hemispheres adjacent to each other with ``view_space_for_other_hemisphere = True``.
      - .. image:: /images/front_template.png
           :width: 400
    * - `medial.h5 <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/view_lookup/medial.h5>`_
      - Correspondence between view from medial side of isocortex and CCF volume. Can place hemispheres near each other with ``view_space_for_other_hemisphere = 'medial'``.
      - .. image:: /images/medial_template.png
           :width: 400
    * - `rotated.h5 <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/view_lookup/rotated.h5>`_
      - Correspondence between view from near the top but rotated to see more of lateral isocortex and CCF volume. Can place hemispheres near each other with ``view_space_for_other_hemisphere = 'rotated'``.
      - .. image:: /images/rotated_template.png
           :width: 400
    * - `side.h5 <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/view_lookup/side.h5>`_
      - Correspondence between view from lateral side of isocortex and CCF volume. Can place hemispheres near each other with ``view_space_for_other_hemisphere = 'side'``.
      - .. image:: /images/side_template.png
           :width: 400
    * - `top.h5 <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/view_lookup/top.h5>`_
      - Correspondence between view from top of isocortex and CCF volume. Can place hemispheres adjacent to each other with ``view_space_for_other_hemisphere = True``.
      - .. image:: /images/top_template.png
           :width: 400


Atlas files (2D)
----------------

.. list-table:: Atlas files
    :widths: 25 50 25
    :header-rows: 1

    * - File
      - Description
      - Image
    * - `back.nrrd <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/master_updated/back.nrrd>`_
      - Atlas of view from back of isocortex
      - .. image:: /images/back_atlas.png
           :width: 400
    * - `bottom.nrrd <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/master_updated/bottom.nrrd>`_
      - Atlas of view from bottom of isocortex
      - .. image:: /images/bottom_atlas.png
           :width: 400
    * - `flatmap_butterfly.nrrd <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/master_updated/flatmap_butterfly.nrrd>`_
      - Atlas of flattened map of all of isocortex - medial side has been adjusted so that the hemispheres are abutting at the center.
      - .. image:: /images/flatmap_butterfly_atlas.png
           :width: 400
    * - `flatmap_dorsal.nrrd <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/master_updated/flatmap_dorsal.nrrd>`_
      - Atlas of flattened map of all of isocortex
      - .. image:: /images/flatmap_dorsal_atlas.png
           :width: 400
    * - `front.nrrd <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/master_updated/front.nrrd>`_
      - Atlas of view from front of isocortex
      - .. image:: /images/front_atlas.png
           :width: 400
    * - `medial.nrrd <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/master_updated/medial.nrrd>`_
      - Atlas of view from medial side of isocortex
      - .. image:: /images/medial_atlas.png
           :width: 400
    * - `rotated.nrrd <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/master_updated/rotated.nrrd>`_
      - Atlas of view from near the top but rotated to see more of lateral isocortex
      - .. image:: /images/rotated_atlas.png
           :width: 400
    * - `side.nrrd <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/master_updated/side.nrrd>`_
      - Atlas of view from lateral side of isocortex
      - .. image:: /images/side_atlas.png
           :width: 400
    * - `top.nrrd <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/master_updated/top.nrrd>`_
      - Atlas of view from top of isocortex
      - .. image:: /images/top_atlas.png
           :width: 400
    * - `labelDescription_ITKSNAPColor.txt <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/master_updated/labelDescription_ITKSNAPColor.txt>`_
      - Metadata for cortical areas in atlas files. Contains identifiers, color, and CCF ontology acronyms.
      -

Isocortex metrics
-----------------

.. list-table:: Metrics files
    :widths: 25 75
    :header-rows: 1

    * - File
      - Description
    * - `avg_layer_depths.json <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/cortical_metrics/avg_layer_depths.json>`_
      - A set of average layer depths to use as targets for normalizing layer thicknesses. **Note:** These were calculated from brain slices in mouse visual cortex - you may want to calculate or use your own thickness values instead.
    * - `cortical_layers_10_v2.h5 <https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/cortical_coordinates/ccf_2017/ccf_streamlines_assets/cortical_metrics/cortical_layers_10_v2.h5>`_
      - The starting depth, ending depth, and thickness of each layer for each streamline. **Warning:** Large file (~100 MB)
