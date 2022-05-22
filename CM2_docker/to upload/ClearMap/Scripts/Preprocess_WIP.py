#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CellMap
=======

This script is the main pipeline to analyze immediate early gene expression 
data from iDISCO+ cleared tissue [Renier2016]_.

See the :ref:`CellMap tutorial </CellMap.ipynb>` for a tutorial and usage.


.. image:: ../Static/cell_abstract_2016.jpg
   :target: https://doi.org/10.1016/j.cell.2020.01.028
   :width: 300

.. figure:: ../Static/CellMap_pipeline.png

  iDISCO+ and ClearMap: A Pipeline for Cell Detection, Registration, and 
  Mapping in Intact Samples Using Light Sheet Microscopy.


References
----------
.. [Renier2016] `Mapping of brain activity by automated volume analysis of immediate early genes. Renier* N, Adams* EL, Kirst* C, Wu* Z, et al. Cell. 2016 165(7):1789-802 <https://doi.org/10.1016/j.cell.2016.05.007>`_
"""
__author__    = 'Christoph Kirst <christoph.kirst.ck@gmail.com>'
__license__   = 'GPLv3 - GNU General Pulic License v3 (see LICENSE)'
__copyright__ = 'Copyright © 2020 by Christoph Kirst'
__webpage__   = 'http://idisco.info'
__download__  = 'http://www.github.com/ChristophKirst/ClearMap2'

#set correct clearmap folder 
import sys
sys.path.append('/home/wirrbel/ClearMap2')

import ClearMap.Environment 

import os, tempfile

if __name__ == "__main__":
     
  #%%############################################################################
  ### Initialization 
  ###############################################################################
  
  #%% Initialize workspace
  
  from ClearMap.Environment import *  #analysis:ignore
  
  #directories and files
  directory = '/home/wirrbel/2020-12-15_CM2_test/Test_imgs_VIP_m3_subset/' #must end on trailing slash   
  
  expression_raw      = '200101_VIP_m_359922_tdtom_19-06-44/14-40-49_VIP_P26_m69368_Het_tdtom_UltraII_C00_xyz-Table Z<Z,4>.ome.tif'           
  expression_auto     = expression_raw
  #expression_auto     = 'Autofluorescence/Auto/Z<Z,4>.tif'  
  
  ws = wsp.Workspace('CellMap', directory=directory);
  ws.update(raw=expression_raw, autofluorescence=expression_auto)
  ws.info()
  
  ws.debug = False
  
  resources_directory = settings.resources_path
  
  #%% Initialize alignment 
  
  #init atlas and reference files
  '''
  annotation_file, reference_file, distance_file=ano.prepare_annotation_files(
      slicing=(slice(None),slice(None),slice(0,256)), 
      orientation=(1,-2,3),
      overwrite=False, verbose=True);
 
 
  #alignment parameter files    
  align_channels_affine_file   = io.join(resources_directory, 'Alignment/align_affine.txt')
  align_reference_affine_file  = io.join(resources_directory, 'Alignment/align_affine.txt')
  align_reference_bspline_file = io.join(resources_directory, 'Alignment/align_bspline.txt')
  '''
  atlas_dir = '/home/wirrbel/2020-12-15_CM2_test/Test_imgs_VIP_m3/atlas_p21'
  atlas_file = io.join(atlas_dir, 'Kim_ref_P21_v2_brain.tif')
  atlas_annotation_file = io.join(atlas_dir,'Kim_ref_P21_v2.1_label.nrrd')
  
  #custom atlas
  annotation_file, reference_file, distance_file = ano.prepare_annotation_files(
      slicing = None,
      directory = atlas_dir,
      annotation_file = atlas_annotation_file,
      reference_file = atlas_file,
      distance_to_surface_file = None,
      orientation = (1,2,3),
      overwrite=False,
      verbose=True)
  
  align_channels_affine_file   = io.join(atlas_dir,'Par0000affine_acquisition.txt')
  align_reference_affine_file  = io.join(atlas_dir,'Par0000affine_acquisition.txt')
  align_reference_bspline_file = io.join(atlas_dir,'Par0000bspline.txt')
  
  #%%############################################################################
  ### Data conversion
  ############################################################################### 
  
  #%% Convet raw data to npy file     
  
  
  #create temp files 
  ws.update(stitched = tempfile.NamedTemporaryFile(prefix = 'stitched', suffix ='.npy').name)  
  ws.update(autofluo = tempfile.NamedTemporaryFile(prefix = 'autofluo', suffix ='.npy').name)
  
  
  source = ws.source('raw');
  sink   = ws.filename('stitched')
  io.delete_file(sink)
  io.convert(source, sink, processes=1, verbose=False, experimentalFixShape=True); #processes = None for all processors
  
  
  #convert autofluo to npy file
  source = ws.source('autofluorescence');
  sink   = ws.filename('autofluo')
  io.delete_file(sink)
  io.convert(source, sink, processes=1, verbose=False, experimentalFixShape=True); #processes = None for all processors
  
  #%%############################################################################
  ### Resampling and atlas alignment 
  ###############################################################################
  
  #%% Resample 
   
  resample_parameter = {
      "source_resolution" : (3.25, 3.25, 3),
      "sink_resolution"   : (25,25,25),
      "processes" : None,
      "verbose" : False,
      "orientation" : (1,2,3)
      };
  
  io.delete_file(ws.filename('resampled'))
  
  res.resample(ws.filename('stitched'), sink=ws.filename('resampled'), **resample_parameter)
  #res.resample(ws.file_list('raw', extension='tif'), sink=ws.filename('resampled'), **resample_parameter)
  
  #%%
  #p3d.plot(ws.filename('resampled'))
  
  #%% Resample autofluorescence
      
  resample_parameter_auto = {
      "source_resolution" : (3.25,3.25,3),
      "sink_resolution"   : (25,25,25),
      "processes" : None,
      "verbose" : False,          
      "orientation": (1,2,3)
      };    
  
  
  res.resample(ws.filename('autofluo'), sink=ws.filename('resampled', postfix='autofluorescence'), **resample_parameter_auto)
  
    #%% remove temp files 
  #because apparently tempfile does not delete them by itself when used with Spyder
  io.remove(ws.filename('stitched'))
  io.remove(ws.filename('autofluo'))