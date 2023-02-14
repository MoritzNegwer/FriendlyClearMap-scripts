# Streamline-related files

/allen/aibs/ccf/2017_integration/handoff

Also copied to /allen/programs/celltypes/workgroups/ivscc/nathang/flatmap/handoff

CCF volume at 10 micron voxels is 1320 x 800 x 1140
* 1320 is x - rostral to caudal
* 800 is y - superior to inferior
* 1140 is z - left to right

## annotation

isocortex_annotation_10.nrrd:
- CCF-sized
- only one hemisphere is labeled
- ids correspond to indices of labels in `master_updated/labelDescription_ITKSNAPColor.txt`

isocortex_annotation_asp.nrrd
- seems slightly smaller and rotated from CCF volume
- **not sure what it's for**


## average_template

back.nrrd
- 1140 x 800 (z by y)
- 2D max projection of average template seen from back
- one hemisphere

bottom.nrrd
- 1140 x 1320 (z by x)
- 2D max projection of average template seen from bottom (assuming empty where streamline end not visible)
- one hemisphere

flatmap_butterfly.nrrd
- 1360 x 1360
- flattened 2D max projection of average template with filled-in notch
- one hemisphere

flatmap_dorsal.nrrd
- 1360 x 1360
- flattened 2D max projection of average template without filled-in notch
- one hemisphere

front.nrrd
- 1140 x 800 (z by y)
- 2D max projection of average template seen from front
- one hemisphere

medial.nrrd
- 1320 x 800 (x by y)
- 2D max projection of average template seen from medial side
- one hemisphere

rotated.nrrd
- 1140 x 1320 (z by x)
- similar to top view but angled to see more of lateral cortex
- one hemisphere

side.nrrd
- 1320 x 800 (x by y)
- 2D max projection of average template seen from lateral side
- one hemisphere

top.nrrd
- 1140 x 1320 (z by x)
- 2D max projection of average template seen from top
- one hemisphere


## cortical_metrics

avg_layer_depths.json

closest_surface_voxel.nrrd

cortical_areas.nrrd

cortical_layers_10.h5
- has keys for each layer
- for each, has one row per path. three values - start distance, end distance, and thickness

cortical_layers.nrrd
- CCF-shaped
- values include 667481440, 667481441, 667481445, 667481446, 667481449, 667481450 (as well as 0)
- these are ontology set IDs (but no way to get their names directly?)
-

depth_from_surface.nrrd

layer_model_depth_mask.nrrd

layer_model_depth.nrrd

layer_relative_depth_mask.nrrd

layer_relative_depth.nrrd

on_path_mask.nrrd

rotation_10.h5

thickness.nrrd


## master

**Presumably older version of `master_updated`??**

filled.nrrd
- **not in `master_updated`**

## master_updated

back.nrrd
- 1140 x 800 (z by y)
- view from back
- ids correspond to indices of labels in `master_updated/labelDescription_ITKSNAPColor.txt`
- one hemisphere

bottom.nrrd
- 1140 x 1320 (z by x)
- view from underneath
- ids correspond to indices of labels in `master_updated/labelDescription_ITKSNAPColor.txt`
- one hemisphere

`flatmap_butterfly.nrrd`
- 1360 x 1360
- flattened isocortex with filled-in notch
- ids correspond to indices of labels in `master_updated/labelDescription_ITKSNAPColor.txt`
- one hemisphere

`flatmap_dorsal.nrrd`
- 1360 x 1360
- flattened isocortex without filled-in notch
- ids correspond to indices of labels in `master_updated/labelDescription_ITKSNAPColor.txt`
- one hemisphere

front.nrrd
- 1140 x 800 (z by y)
- view from front
- ids correspond to indices of labels in `master_updated/labelDescription_ITKSNAPColor.txt`
- one hemisphere

labelDescription_ITKSNAPColor.txt
- list of indices, colors, and acronyms of structures

master.nrrd
- 1380 x 800 x 1140 (CCF shaped)
- ids correspond to indices of labels in `master_updated/labelDescription_ITKSNAPColor.txt`
- only tops of streamlines are annotated with structures
- can be used to create different 2D views
- one hemisphere

medial.nrrd
- 1320 x 800 (x by y)
- view from medial side
- ids correspond to indices of labels in `master_updated/labelDescription_ITKSNAPColor.txt`
- one hemisphere

rotated.nrrd
- 1140 x 1320 (z by x)
- similar to top view but angled to see more of lateral cortex
- ids correspond to indices of labels in `master_updated/labelDescription_ITKSNAPColor.txt`
- one hemisphere

side.nrrd
- 1320 x 800 (x by y)
- view from lateral side
- ids correspond to indices of labels in `master_updated/labelDescription_ITKSNAPColor.txt`
- one hemisphere

top.nrrd
- 1140 x 1320 (z by x)
- view from top
- ids correspond to indices of labels in `master_updated/labelDescription_ITKSNAPColor.txt`
- one hemisphere



## streamlines

surface_paths_10.h5
- key: "paths"
- shape (1476024, 200)
- contains flattened voxel locations of members of path (up to 200 elements; unused ones are set to zero)

- key: "volume lookup"
- shape (1320, 800, 1140) (CCF shape)
- volume has the index of the path for the voxel

isocortex_boundary_10.nrrd
- CCF-shaped
- contains values from 0 to 4 **what do they represent? - maybe top, bottom, etc?**
- one hemisphere

isocortex_layers_10.nrrd
- CCF-shaped
- one hemisphere
- filled in isocortex with indices indicated layer (uses the "VIS" set in the labels file for the entire cortex)

isocortex_mask_10.nrrd
- CCF-shaped
- one hemisphere
- filled in isocortex mask - value of 2 means inside isocortex

isocortex_top_surface_10.nrrd
- CCF-shaped
- one hemisphere
- values of 1 for voxels belonging to top of cortex

closest_surface_voxel_index.nrrd
- CCF-shaped
- one hemisphere
- contains flattened index of best-match surface voxel
- useful for finding best streamline for voxel

laplacian_10.nrrd
- CCF-shaped
- one hemisphere
- value of the laplacian between top and bottom of isocortex


## view_lookup

back.h5
- contains dataset "view lookup"
- has attrs 'origin', 'size' (CCF size), 'spacing' (x-y-z resolution in microns), 'view size' - size of output
- first column is flattened indices for the back 2D projection
- second column is the flattened indices of the CCF volume for the 'volume lookup' data structure
- presumably second column is like this is because this is easy to generate? it's just connecting the 2D location to the point on the 3D volume that you are viewing


# Basic operations

to do max projection:
1. take the projection h5 file (e.g. "back.h5")
2. use second column of 'view lookup' with the surface_paths_10.h5 'volume lookup' to find the indices of the paths to use
3. grab just those paths from 'paths' in surface_paths_10.h5 (and order them to match 'view lookup')
4. for each path, grab voxels from flattened input and get maximum (using correct row of paths) and assign to correct flattened index of output

