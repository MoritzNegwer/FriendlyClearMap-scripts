#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 23:47:22 2022

@author: wirrbel
"""



import nrrd 
import os
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
import tifffile

import sys
sys.path.append('../cortical_flatmaps/ccf_streamlines-main/')
import ccf_streamlines.projection as ccfproj

data_dir = '../cortical_flatmaps/'
output_dir = '../cortical_flatmaps/flatmaps_test/'

#if output_dir doesn't exist, create it
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def parseElastixOutputPoints(filename, indices = False):
    #from Clearmap 1 Elastix.py
    """Parses the output points from the output file of transformix

    Arguments:
        filename (str): file name of the transformix output file
        indices (bool): if True return pixel indices otherwise float coordinates

    Returns:
        points (array): the transformed coordinates     
    """

    with open(filename) as f:
        lines = f.readlines()
        f.close();

    length = len(lines);

    if length == 0:
        return np.zeros((0,3));

    points = np.zeros((length, 3));
    k = 0;
    for line in lines:
        ls = line.split();
        if indices:
            for i in range(0,3):
                points[k,i] = float(ls[i+22]);
        else:
            for i in range(0,3):
                points[k,i] = float(ls[i+30]);

        k += 1;

    return points;


def summarize_all_cells_per_group_as_10_um_nrrd(mouse_list,output_dir,output_name):
    all_cells = np.empty(shape=(0,3),dtype="float32")
    
    for mouse in mouse_list:
        cells_dir = mouse[0]
        #load CCF3-transformed cells as floating-point px coordinates 
        cells = parseElastixOutputPoints(os.path.join(cells_dir,'Aligned_CCF3/outputpoints.txt'),indices=False).astype('float32')
        #add to all_cells
        all_cells = np.append(all_cells, cells,axis=0)
    

    
    
    #resample to 10um spacing (avg 10 um atlas has coordinates x=1320, y=800, z=1140)
    #25 um atlas x = 528, y = 320 , z = 456
    all_cells = (np.array(all_cells)*2.5).astype('float32')
    
    #define heatmap with xyz-coords of 10um CCF atlas
    heatmap = np.zeros(shape=(1320,800,1140)).astype('uint16')
    
    #add each cell to heatmap 
    skipped_cells = 0
    for cell in all_cells:
        try:
            heatmap[round(cell[0]),round(cell[1]),round(cell[2])] += 1
        except: 
            skipped_cells += 1
            pass
    print ('cells skipped for being out of bounds: ',skipped_cells)
    
    #normalize for # of samples in group
    heatmap = heatmap.astype('float32')/len(mouse_list)
    
    #blur heatmap 
    heatmap = gaussian_filter(heatmap, sigma=5,output='float32')
        
    #TODO: save as nrrd 
    output_file_path = os.path.join(output_dir,output_name+'.nrrd')
    nrrd.write(output_file_path,heatmap)
    
    #debug: save as tiff stack
    #save tiff stack
    tifffile.imwrite(os.path.join(output_dir,output_name+'.tif'),heatmap,compression="zlib")
    
    return output_file_path
    
    

#==== P14 ====
WorkFoldersP14 = [['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P14/2020-03-02_Lynn_PV_P14_m1/','WT',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P14/2020-03-02_Lynn_PV_P14_m2/','WT',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P14/2020-03-02_Lynn_PV_P14_m3/','WT',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P14/2020-03-02_Lynn_PV_P14_m5/','WT',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lukas_PV_P14/Lukas_PV_P14_m1_66632/','WT',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lukas_PV_P14/Lukas_PV_P14_m2_66633/','WT',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lukas_PV_P14/Lukas_PV_P14_m6_190926/','WT',14]]

PV_P14_nrrd = summarize_all_cells_per_group_as_10_um_nrrd(WorkFoldersP14,output_dir,'PV_P14_WT')

  
#==== P56 ==== 
WorkFoldersP56 = [['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/Lukas_PV_50_m6/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/Lukas_PV_50_m7/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/Lukas_PV_50_m8/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/2020-02-17_Lynn_PV_P50_m1/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m4/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m6/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/2020-02-24_Lynn_PV_P50_m7/','WT',56],
               ]
    
PV_P56_nrrd = summarize_all_cells_per_group_as_10_um_nrrd(WorkFoldersP56,output_dir,'PV_P56_WT')

#===SST===
WorkFoldersSST = [["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/SST_Maren/2019-04-24_SST_P25_Sample2_57646_WT/", 'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/SST_Maren/2019-04-25_SST_P25_Sample3_57647_WT/", 'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/SST_Maren/2019-04-26_SST_P25_Sample5_57654_WT/", 'WT',28], #WT
               ]
    
SST_P28_nrrd = summarize_all_cells_per_group_as_10_um_nrrd(WorkFoldersSST,output_dir,'SST_P28_WT')

#===VIP===
WorkFoldersVIP = [["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_docker/Lynn_VIP_P26_m69366_WT/", 'WT', 28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_docker/Lynn_VIP_P26_m69369_WT/VIP_P26_m69369_WT_tdtom/",'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_docker/Lynn_VIP_P26_m69370_WT/",'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_docker/Lynn_VIP_P26_m69371_WT/", 'WT', 28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_Arivis/m1/", 'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_Arivis/m3/", 'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_Arivis/m5/", 'WT',28], #WT
                ]

VIP_P28_nnrd = summarize_all_cells_per_group_as_10_um_nrrd(WorkFoldersVIP,output_dir,'VIP_P28_WT')




#### FLATMAP #### 
#adapted from the example code: https://ccf-streamlines.readthedocs.io/en/latest/guide.html

#load template 
template, _ = nrrd.read(os.path.join(data_dir,"average_template_10.nrrd"))

#make butterfly projection 
proj_bf = ccfproj.Isocortex2dProjector(
    # Specify our view lookup file
    os.path.join(data_dir,"flatmap_butterfly.h5"),

    # Specify our streamline file
    os.path.join(data_dir,"surface_paths_10_v3.h5"),

    # Specify that we want to project both hemispheres
    hemisphere="both",

    # The butterfly view doesn't contain space for the right hemisphere,
    # but the projector knows where to put the right hemisphere data so
    # the two hemispheres are adjacent if we specify that we're using the
    # butterfly flatmap
    view_space_for_other_hemisphere='flatmap_butterfly',
)

bf_projection_max = proj_bf.project_volume(template)

'''
plt.imshow(
    bf_projection_max.T, # transpose so that the rostral/caudal direction is up/down
    interpolation='none',
    cmap='Greys_r',
)
'''

#generate boundaries 
bf_boundary_finder = ccfproj.BoundaryFinder(
    projected_atlas_file = os.path.join(data_dir,"flatmap_butterfly.nrrd"),
    labels_file = os.path.join(data_dir,"labelDescription_ITKSNAPColor.txt"),
)

# We get the left hemisphere region boundaries with the default arguments
bf_left_boundaries = bf_boundary_finder.region_boundaries()

# And we can get the right hemisphere boundaries that match up with
# our projection if we specify the same configuration
bf_right_boundaries = bf_boundary_finder.region_boundaries(
    # we want the right hemisphere boundaries, but located in the right place
    # to plot both hemispheres at the same time
    hemisphere='right_for_both',

    # we also want the hemispheres to be adjacent
    view_space_for_other_hemisphere='flatmap_butterfly',
)
#optional: plot the flatmap
'''
plt.imshow(
    bf_projection_max.T,
    interpolation='none',
    cmap='Greys_r',
)

for k, boundary_coords in bf_left_boundaries.items():
    plt.plot(*boundary_coords.T, c="white", lw=0.5)
for k, boundary_coords in bf_right_boundaries.items():
    plt.plot(*boundary_coords.T, c="white", lw=0.5)
''' 

def calculate_cortical_depths (tlx_normalized_layers,layer_tops):
    #added: calculate the normalized_layer depth in voxels 
    depth_voxels = tlx_normalized_layers.shape[2]
    max_cortex_depth_um = layer_tops['wm']
    um_to_voxels = max_cortex_depth_um/depth_voxels
    
    top_l1 = 0
    top_l23 = int(layer_tops['2/3']/um_to_voxels)
    top_l4 = int(layer_tops['4']/um_to_voxels)
    top_l5 = int(layer_tops['5']/um_to_voxels)
    top_l6 = int(layer_tops['6a']/um_to_voxels)
    top_wm = int(layer_tops['wm']/um_to_voxels)
    
    
    #create a list of all layers you want to project 
    layers_flatmap_list = [
    [top_l1,top_l23,"Layer 1"],
    [top_l23,top_l4,"Layer 2/3"],
    [top_l4,top_l5,"Layer 4"],
    [top_l5,top_l6,"Layer 5"],
    [top_l6,top_wm,"Layer 6"],
    [top_l1,top_wm,"Cortex"]]

    return layers_flatmap_list

#### Make Heatmap projection ####
def make_heatmap(data_dir,output_dir,input_data):
    #load density data 
    tlx_data, _ = nrrd.read(input_data)
    
    ### Layer-normalized projection
    #configure projector
    proj_butterfly_slab = ccfproj.Isocortex3dProjector(
        # Similar inputs as the 2d version...
        os.path.join(data_dir,"flatmap_butterfly.h5"),
        os.path.join(data_dir,"surface_paths_10_v3.h5"),
    
        hemisphere="both",
        view_space_for_other_hemisphere='flatmap_butterfly',
    
        # Additional information for thickness calculations
        thickness_type="normalized_layers", # each layer will have the same thickness everwhere
        layer_thicknesses=layer_thicknesses,
        streamline_layer_thickness_file=os.path.join(data_dir,"cortical_layers_10_v2.h5"),
    )
    
    #project data
    #n.b. "normalized layers" = all layers are equally thick, useful to get layer-separated plots 
    tlx_normalized_layers = proj_butterfly_slab.project_volume(tlx_data)
    
    
    ### unnormalized view data for depth profile + tiff export
    
    #configure projector
    proj_butterfly_slab_unnormalized = ccfproj.Isocortex3dProjector(
        # Similar inputs as the 2d version...
        os.path.join(data_dir,"flatmap_butterfly.h5"),
        os.path.join(data_dir,"surface_paths_10_v3.h5"),
    
        hemisphere="both",
        view_space_for_other_hemisphere='flatmap_butterfly',
    
        # Additional information for thickness calculations
        thickness_type="unnormalized", # thickness 
        layer_thicknesses=layer_thicknesses,
        streamline_layer_thickness_file=os.path.join(data_dir,"cortical_layers_10_v2.h5"),
    )
    
    
    tlx_unnormalized = proj_butterfly_slab_unnormalized.project_volume(tlx_data)
    
    #optional: plot cortical layers as grid 
    '''
    # Calculate our maximum intensity projections of the slab
    # and keep track of their shapes (to set up the plots to be right size
    
    main_max = tlx_normalized_layers.max(axis=2).T
    top_max = tlx_normalized_layers.max(axis=1).T
    left_max = tlx_normalized_layers.max(axis=0)
    
    main_shape = main_max.shape
    top_shape = top_max.shape
    left_shape = left_max.shape
    
    
    # Set up a figure to plot them together
    fig, axes = plt.subplots(2, 2,
                             gridspec_kw=dict(
                                 width_ratios=(left_shape[1], main_shape[1]),
                                 height_ratios=(top_shape[0], main_shape[0]),
                                 hspace=0.01,
                                 wspace=0.01),
                             figsize=(19.4, 12))
    
    # Plot the surface view
    axes[1, 1].imshow(main_max, vmin=0, vmax=1, cmap="magma", interpolation=None)
    
    # plot our region boundaries
    for k, boundary_coords in bf_left_boundaries.items():
        axes[1, 1].plot(*boundary_coords.T, c="white", lw=0.5)
    for k, boundary_coords in bf_right_boundaries.items():
        axes[1, 1].plot(*boundary_coords.T, c="white", lw=0.5)
    
    axes[1, 1].set(xticks=[], yticks=[], anchor="NW")
    
    # Plot the top view
    axes[0, 1].imshow(top_max, vmin=0, vmax=1, cmap="magma", interpolation=None)
    axes[0, 1].set(xticks=[], yticks=[], anchor="SW")
    
    # Plot the side view
    axes[1, 0].imshow(left_max, vmin=0, vmax=1, cmap="magma", interpolation=None)
    axes[1, 0].set(xticks=[], yticks=[], anchor="NE")
    
    # Remove axes from unused plot area
    axes[0, 0].set(xticks=[], yticks=[])
    '''
    
    #calculate layer start list 
    layers_flatmap_list = calculate_cortical_depths(tlx_normalized_layers, layer_tops)
    
    return tlx_normalized_layers, layers_flatmap_list, tlx_unnormalized




def layer_flatmap(tlx_normalized_layers, bf_left_boundaries, bf_right_boundaries, layer_top,layer_bottom,heatmap_name,layer_title,output_dir,scale_colormap=True,colormap="inferno"):
    
    #scale the colormap to the data, otherwise, use 1 (Note that the heatmap input might contain values much higher than that)
    if scale_colormap:
        vmax = tlx_normalized_layers.max()
    else:
        vmax = 0.01
    
    # find the max projection of just one layer 
    plt.figure(figsize=(10,10),dpi=300)
    plt.axis('off')
    plt.imshow(
        tlx_normalized_layers[:, :, layer_top:layer_bottom].max(axis=2).T,
        vmin=0, vmax=vmax,
        cmap=colormap,
        interpolation=None
    )
    
    # plot region boundaries
    for k, boundary_coords in bf_left_boundaries.items():
        plt.plot(*boundary_coords.T, c="white", lw=0.1)
    for k, boundary_coords in bf_right_boundaries.items():
        plt.plot(*boundary_coords.T, c="white", lw=0.1)
    
    plt.title(layer_title)
    #save figure as 2D PNG    
    plt.savefig(os.path.join(output_dir,heatmap_name+"_Flatmap_"+layer_title.replace('/','')+".png"))

def save_flatmap_stack(tlx_normalized_layers,output_dir,heatmap_name):
    #flip the heatmap so Z is cortical depth
    flat_viewed_depthwise = tlx_normalized_layers.T
    #save as tif 
    tifffile.imwrite(os.path.join(output_dir,heatmap_name+"_flatmap.tif"),flat_viewed_depthwise,compression="lzw")
    
def generate_depth_profiles (input_collection,output_dir,max_value=False):
    #plots cortical depth profiles
    
    #create collection 
    depth_profile_collection = {}
    
    #fill the collection with depth maps, save individual lists
    for heatmap_name, heatmap in input_collection.items():
        #reduce to 1 dimension (=cortical depth)
        depth_profile = np.sum(heatmap,axis=0)
        depth_profile = np.sum(depth_profile,axis=0)
      
        #add to collection
        depth_profile_collection[heatmap_name] = depth_profile
      
        
        #create single plot, then save 
        plt.clf()
        plt.plot(depth_profile)
        plt.title(heatmap_name)
        plt.savefig(os.path.join(output_dir,"depth_map_"+heatmap_name+"_v01.svg"))  
        
    #create overview with all on the same scaling
    fig, axes = plt.subplots(len(depth_profile_collection), 1,figsize=(10,10),sharex=True)
    #set up to display depth in µm (rather than count of 10µm voxels)
    plt.setp(axes, xticks=[0,100,200],xticklabels=[0,1000,2000])
    #add labels
    fig.supxlabel("Cortical Depth (µm)")
    fig.supylabel("Normalized Cell Count (n)")
    #if a max value is supplied, use it here 
    if max_value is not False:
        plt.setp(axes,yticks=[0,int(round(max_value/2)),max_value],ylim=[0,max_value])
    
    #plot all cells 
    for i, heatmap_name in enumerate(depth_profile_collection):
    
        axes[i].plot(depth_profile_collection[heatmap_name])
        #axes[i].set_title(heatmap_name)
        axes[i].set_ylabel(heatmap_name)
        
    #save as svg
    fig.savefig(os.path.join(output_dir,"depth_map_collection_all_v01.svg"))
    
    
#### PREPARATION 

#calculate layer depths 
# Note - the layer depths in this file was calculated from a set of visual cortex slices,
# so you may wish to use another set of layer depths depending on your purposes
with open(os.path.join(data_dir,"avg_layer_depths.json"), "r") as f:
    layer_tops = json.load(f)

layer_thicknesses = {
        'Isocortex layer 1': layer_tops['2/3'],
        'Isocortex layer 2/3': layer_tops['4'] - layer_tops['2/3'],
        'Isocortex layer 4': layer_tops['5'] - layer_tops['4'],
        'Isocortex layer 5': layer_tops['6a'] - layer_tops['5'],
        'Isocortex layer 6a': layer_tops['6b'] - layer_tops['6a'],
        'Isocortex layer 6b': layer_tops['wm'] - layer_tops['6b'],
}

### Plot single layers 

all_nrrds = [[PV_P14_nrrd,'PV_P14'],
             [PV_P56_nrrd,'PV_P56'],
             [SST_P28_nrrd,'SST_P28'],
             [VIP_P28_nnrd,'VIP_P28'],
             ]

#for depth map: generate depth map collection
depth_map_collection = {}

#go through all nrrd files and create tlx_normalized_layers 
for nrrd_file,heatmap_name in all_nrrds:
    #create tlx_normalized_layers 
    tlx_normalized_layers, layers_flatmap_list, tlx_unnormalized = make_heatmap(data_dir,output_dir,nrrd_file)
    
    #optional:save tlx_normalized_layers as tiff (true flatmap)
    save_flatmap_stack(tlx_normalized_layers,output_dir,heatmap_name+"_equal_layer_thickness")
    
    #optional: make depth histogram + save tiff
    save_flatmap_stack(tlx_unnormalized,output_dir,heatmap_name+"_original_cortex_thickness")
    depth_map_collection[heatmap_name] = tlx_unnormalized
    
    #save one plot each 
    for layer_top,layer_bottom,layer_title in layers_flatmap_list:
        layer_flatmap(tlx_normalized_layers, bf_left_boundaries, bf_right_boundaries, layer_top,layer_bottom,heatmap_name,layer_title,output_dir,scale_colormap=False,colormap="inferno")
    
#compute all depth profiles and save 
generate_depth_profiles(depth_map_collection, output_dir,max_value=8000)

