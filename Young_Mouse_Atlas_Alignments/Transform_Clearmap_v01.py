#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 16:09:13 2021

@author: wirrbel
"""

####Testing alignment with ClearMap components 
import sys
sys.path.append('/home/wirrbel/ClearMap2')

from ClearMap.Environment import *  #analysis:ignore

#directories and files
directory = '/home/wirrbel/2021-03-29_brainrender_2_preprocessing/CCF3_Clearmap/' #must end on trailing slash

ws = wsp.Workspace('CellMap', directory=directory);

ws.debug = False

resources_directory = settings.resources_path

atlas_dir = '/home/wirrbel/2021-03-29_brainrender_2_preprocessing/Kim_dev_atlas/'
atlas_file = io.join(atlas_dir, 'CCF3_average_template_25.tif')
moving_file = io.join(atlas_dir, 'Kim_ref_P14_v2_brain.tif')

annotation_file, reference_file, distance_file = ano.prepare_annotation_files(
    slicing = None,
    directory = atlas_dir,
    annotation_file = atlas_file,
    reference_file = atlas_file,
    distance_to_surface_file = None,
    orientation = (1,2,3),
    overwrite=False,
    verbose=True)

align_channels_affine_file   = io.join('/home/wirrbel/2021-03-29_brainrender_2_preprocessing/elastix_files/','Par0000affine.txt')
align_reference_affine_file  = io.join('/home/wirrbel/2021-03-29_brainrender_2_preprocessing/elastix_files/','Par0000affine.txt')
align_reference_bspline_file = io.join('/home/wirrbel/2021-03-29_brainrender_2_preprocessing/elastix_files/','Par0000bspline.txt')
  

# align autofluorescence to reference
align_reference_parameter = {            
    #moving and reference images
    "moving_image" : reference_file,
    "fixed_image"  : moving_file,
    
    #elastix parameter files for alignment
    "affine_parameter_file"  :  align_reference_affine_file, #Note that this pre-alignment is mandatory for now (2020-12-15)
    "bspline_parameter_file" :  align_reference_bspline_file,
    #directory of the alignment result
    "result_directory" :  ws.filename('auto_to_reference'),
    #for BigWarp manual alignment 
    "moving_points" : io.join(directory,'atlas_landmarks.txt'), #added by Moritz 2020-10-20
    "fixed_points" : io.join(directory,'autofluo_landmarks.txt') #added by Moritz 2020-10-20
    };

elx.align(**align_reference_parameter);


### === load cells 
# Resample autofluorescence
'''
resample_parameter_auto = {
    "source_resolution" : (1,1,1),
    "sink_resolution"   : (1,1,1),
    "processes" : None,
    "verbose" : False,          
    "orientation": (1,2,3)
    };    


res.resample(ws.filename('autofluo'), sink=ws.filename('resampled', postfix='autofluorescence'), **resample_parameter_auto)
'''



####==== Transform cells 

#source = ws.source('cells', postfix='filtered')
source = np.load('/home/wirrbel/2021-03-29_brainrender_2_preprocessing/CCF3_Clearmap/cells_transformed_to_Atlas_aligned.npy')

def transformation(coordinates):
  #debug
  print ("starting resampling points, coords: ",coordinates[0])
  '''
  coordinates = res.resample_points(
                  coordinates, sink=None, orientation=None, 
                  source_shape=io.shape(ws.filename('stitched')), 
                  sink_shape=io.shape(ws.filename('resampled')),
   
                  );
  
  #debug
  print ("starting transforming points resamnpled -> auto, coords: ",coordinates[0])
  
  coordinates = elx.transform_points(
                  coordinates, sink=None, 
                  transform_directory=ws.filename('resampled_to_auto'), 
                  binary=False, indices=False);
  
  #debug
  print ("starting transforming points auto -> reference, coords: ",coordinates[0])
  '''
  coordinates = elx.transform_points(
                  coordinates, sink=None, 
                  transform_directory=ws.filename('auto_to_reference'),
                  binary=False, indices=False);
      
  return coordinates;
  

coordinates = source;

coordinates_transformed = transformation(coordinates);

