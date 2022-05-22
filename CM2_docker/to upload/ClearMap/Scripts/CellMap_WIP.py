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
  directory = '/home/wirrbel/2020-12-15_CM2_test/Test_imgs_VIP_m3/' #must end on trailing slash   
  
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
  

  #p3d.plot([ws.filename('resampled'), ws.filename('resampled', postfix='autofluorescence')])
  
  #%% Alignment - resampled to autofluorescence
  
  # align the two channels
  align_channels_parameter = {            
      #moving and reference images
      "moving_image" : ws.filename('resampled', postfix='autofluorescence'),
      "fixed_image"  : ws.filename('resampled'),
      
      #elastix parameter files for alignment
      "affine_parameter_file"  : align_channels_affine_file,
      "bspline_parameter_file" : None,
      
      #directory of the alig'/home/nicolas.renier/Documents/ClearMap_Ressources/Par0000affine.txt',nment result
      "result_directory" :  ws.filename('resampled_to_auto')
      
      }; 
  
  elx.align(**align_channels_parameter);
  
  #%% Alignment - autoflourescence to reference
  
  # align autofluorescence to reference
  align_reference_parameter = {            
      #moving and reference images
      "moving_image" : reference_file,
      "fixed_image"  : ws.filename('resampled', postfix='autofluorescence'),
      
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
  
  
  #%%############################################################################
  ### Create test data
  ###############################################################################
  
  #%% Crop test data 
  
  #select sublice for testing the pipeline
  #slicing = (slice(2000,2200),slice(2000,2200),slice(50,150));
  #ws.create_debug('stitched', slicing=slicing);
  #ws.debug = True; 
  
  #p3d.plot(ws.filename('stitched'))
    
  
  #%%############################################################################
  ### Cell detection
  ###############################################################################
  
  #%% Cell detection:
  
  cell_detection_parameter = cells.default_cell_detection_parameter.copy();
  cell_detection_parameter['illumination_correction'] = None;
  cell_detection_parameter['background_correction'] = None;
  cell_detection_parameter['intensity_detection']['measure'] = ['source'];
  cell_detection_parameter['shape_detection']['threshold'] = 500;
  
  #io.delete_file(ws.filename('cells', postfix='maxima'))
  #cell_detection_parameter['maxima_detection']['save'] = ws.filename('cells', postfix='maxima')
  
  processing_parameter = cells.default_cell_detection_processing_parameter.copy();
  processing_parameter.update(
      processes = None, # Ilastik
      size_max = 300, #100, #35,
      size_min = 20,# 30, #30,
      overlap  = 5, #32, #10,
      verbose = True
      )
  
  ##Optional: Cell detection with Ilastik ported over from CM1####
  #old:
  #ImageProcessingMethod = "Ilastik";
  cell_detection_parameter['cell_detection_method'] = 'Ilastik'
  
  #add IlastikParameter and findCellsIntensity Parameter to cell_detection_parameter
  cell_detection_parameter['classifyCellsParameter'] = dict(
      classifier = "/home/wirrbel/2020-11-26_WFA_Ilastik/2020-11-26_WFA_Ilastik_03_8bitexport.ilp",
      classindex = 0,
      save      = None, #os.path.join('/home/wirrbel/2020-11-26_WFA_Ilastik/temp_ilastik_save/','IlastikClasses/Z\d{4}.ome.tif'), #os.path.join(BaseDirectory, 'IlastikClasses/Z\d{4}.ome.tif'),      # (str or None)       file name to save result of this operation if None dont save to file 
      processingDirectory  = None,  #where the temp file will be placed. If this is set to None, the file will be generated in /tmp, i.e. the RAM. If you have <256GB RAM, it will likely fill the memory and trigger the OOM-killer. 
      verbose = True       # (bool or int)       print / plot information about this step
      )
  
  cell_detection_parameter['findCellIntensityParameter'] = dict(
      method = 'Max',  # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
      verbose = True, # (bool or int)       print / plot information about this step
      intensities_sink = io.join(directory,'intensities_raw.npy')
      )
  
  cells.detect_cells(source = ws.source('raw'), 
                     sink = ws.filename('cells', postfix='raw'),
                     cell_detection_parameter=cell_detection_parameter, 
                     processing_parameter=processing_parameter)
  

  
  
  
  #%% visualization
  
  #p3d.plot([[ws.filename('stitched'), ws.filename('cells', postfix='maxima')]])
  
  #%%
  #coordinates = np.hstack([ws.source('cells', postfix='raw')[c][:,None] for c in 'xyz']);
  #p = p3d.list_plot_3d(coordinates)
  #p3d.plot_3d(ws.filename('stitched'), view=p, cmap=p3d.grays_alpha(alpha=1))
  
  
  #%% Cell statistics
  '''
  source = ws.source('cells', postfix='raw')
  #source = io.join(directory,'cells_allpoints_aligned.npy')
  
  plt.figure(1); plt.clf();
  names = source.dtype.names;
  nx,ny = p3d.subplot_tiling(len(names));
  for i, name in enumerate(names):
    plt.subplot(nx, ny, i+1)
    plt.hist(source[name]);
    plt.title(name)
  plt.tight_layout();
  '''
  #%% Filter cells
  
  
  thresholds = {
      'source' : (200,None), #original: 'source': None
      'size'   : None
      }
  
  cells.filter_cells(source = ws.filename('cells', postfix='raw'), 
                     sink = ws.filename('cells', postfix='filtered'), 
                     thresholds=thresholds);
  
  #cells.filter_cells(source = io.join(directory,'cells_allpoints_aligned.npy'),
  #                    sink = ws.filename('cells', postfix='filtered'), 
  #                    thresholds=thresholds);
  
  
  #%% Visualize
  
  #coordinates = np.array([ws.source('cells', postfix='filtered')[c] for c in 'xyz']).T;
  #p = p3d.list_plot_3d(coordinates, color=(1,0,0,0.5), size=10)
  #p3d.plot_3d(ws.filename('stitched'), view=p, cmap=p3d.grays_alpha(alpha=1))
  
  
  #%%############################################################################
  ### Cell atlas alignment and annotation
  ###############################################################################
  
  #%% Cell alignment
  
  #source = ws.source('cells', postfix='filtered')
  source = ws.source('cells', postfix='filtered')
  
  def transformation(coordinates):
    #debug
    print ("starting resampling points, coords: ",coordinates[0])
    
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
    
    coordinates = elx.transform_points(
                    coordinates, sink=None, 
                    transform_directory=ws.filename('auto_to_reference'),
                    binary=False, indices=False);
        
    return coordinates;
    
  
  coordinates = np.array([source[c] for c in 'xyz']).T;
  
  coordinates_transformed = transformation(coordinates);
  
  #debug
  print ("transformed coords: ",coordinates_transformed)
  
  #%% Cell annotation
  #load custom label file 
  import pandas
  custom_label_csv = pandas.read_csv('/home/wirrbel/2020-12-15_CM2_test/Test_imgs_VIP_m3/atlas_p21/regions_IDs.csv')
  #load default extra label 
  extra_label = ano.default_extra_label
  
  def add_missing_value (missing_value,extra_label):
    #find missing value in custom file 
    value_to_add = custom_label_csv.loc[custom_label_csv['id'] == missing_value]
    
    if value_to_add.empty: 
        #if the value is not in the file, add as background 
        value_to_add = custom_label_csv.loc[custom_label_csv['id'] == 0]
        value_to_add['parent_id'] = 997
    else:
        #if value is in custom_label_csv, add 
        #but first check if parent is present 
        parent_value = int(value_to_add['parent_id'].item())
        try:
            ano.find(parent_value,key='id')
        except KeyError:
            #if the parent id isn't present, add recursively
            print ("addidg missing parent value: ",parent_value)
            extra_label = add_missing_value(parent_value,extra_label)
                
    #format into extra label
    extra_label.append((missing_value,int(value_to_add['parent_id'].item()),value_to_add['name'].item(),value_to_add['acronym'].item()))
    
    print ("value added: ", (missing_value,int(value_to_add['parent_id'].item()),value_to_add['name'].item(),value_to_add['acronym'].item()))

    return extra_label
  
  
  #added annotation file to make sure that custom atlases are taken along 
  def try_label(coordinates_transformed = coordinates_transformed,annotation_file = annotation_file, extra_label = extra_label):
      try: 
          print ("trying to label all cells")
          label = ano.label_points(coordinates_transformed, annotation_file = annotation_file, key='order');
          names = ano.convert_label(label, key='order', value='name');
          success = True
      except KeyError as missing_value:
          #make number
          missing_value = missing_value.args[0]         
          print ("Key error with key: ", missing_value)  
          #recursively add missing_value and its parents (if not present)
          extra_label = add_missing_value(missing_value,extra_label)
          
          success = False
          
          #re-initialize with new value
          
          #run again, hopefully with extra label
          
          #success, label, extra_label = try_label(coordinates_transformed, annotation_file = annotation_file, extra_label = extra_label)
      
      ano.set_label_file(label_file = ano.default_label_file, extra_label=extra_label)
      if not success:
          success, label, extra_label, names = try_label(coordinates_transformed, annotation_file = annotation_file, extra_label = extra_label)
      #return label, extra_label
      
      return success, label, extra_label, names
   
  success = False
  
  while not success:
      success, label, extra_label,names = try_label (coordinates_transformed = coordinates_transformed,annotation_file = annotation_file, extra_label = extra_label)
  #ano.set_label_file(label_file = ano.default_label_file, extra_label=extra_label)
  #label = ano.label_points(coordinates_transformed, annotation_file = annotation_file, key='order');
  
  #debug
  print("all labels appended, extracting names")
  
  
  #names = ano.convert_label(label, key='id', value='name');
  
  #%% Save results
  
  coordinates_transformed.dtype=[(t,float) for t in ('xt','yt','zt')]
  label = np.array(label, dtype=[('order', int)]);
  names = np.array(names, dtype=[('name' , 'U256')])
  
  import numpy.lib.recfunctions as rfn
  cells_data = rfn.merge_arrays([source[:], coordinates_transformed, label, names], flatten=True, usemask=False)
  
  io.write(ws.filename('cells'), cells_data)
  
  
  
  #%%############################################################################
  ### Cell csv generation for external analysis
  ###############################################################################
  
  #%% CSV export
  
  source = ws.source('cells');
  header = ', '.join([h[0] for h in source.dtype.names]);
  np.savetxt(ws.filename('cells', extension='csv'), source[:], header=header, delimiter=',', fmt='%s')
  
  #%% ClearMap 1.0 export
  
  source = ws.source('cells');
  
  clearmap1_format = {'points' : ['x', 'y', 'z'], 
                      'points_transformed' : ['xt', 'yt', 'zt'],
                      'intensities' : ['source', 'dog', 'background', 'size']}
  
  for filename, names in clearmap1_format.items():
    sink = ws.filename('cells', postfix=['ClearMap1', filename]);
    data = np.array([source[name] if name in source.dtype.names else np.full(source.shape[0], np.nan) for name in names]);
    io.write(sink, data);
    
  
  #CM1-style table export
  #TODO: add extra label 
  counts = ano.count_label(label['order'],weights=None,key='order', hierarchical=True)
  ids_list = ano.get_list('id')
  name_list = ano.get_list('name')

  #additional info 
  acronym_list = ano.get_list('acronym')
  order_list = ano.get_list('graph_order')
  parent_list = ano.get_list('parent_structure_id')
  level_list = ano.get_list('level')
  color_list = ano.get_list('color_hex_triplet')


  table = np.zeros(counts.shape, dtype=[('id','int64'),('counts','f8'),('name', 'U256'),('acronym', 'U256'),('order','U256'),('parent_structure_id','U256'),('level','U256'),('color_hex_triplet','U256')])
  table['id'] = ids_list
  table['counts'] = counts
  table['name'] = name_list

  #additional info 
  table['acronym'] = acronym_list
  table['order'] = order_list
  table['parent_structure_id'] = parent_list
  table['level'] = level_list
  table['color_hex_triplet'] = color_list

  #export to csv file 
  np.savetxt(io.join(directory, 'Annotated_counts_ClearMap_1.csv'), table, delimiter=';', fmt='%s',
             header="id; counts;name;acronym;order;parent_structure_id;level;color_hex_triplet")  

    
  #%%############################################################################
  ### Voxelization - cell density
  ###############################################################################
  
  source = ws.source('cells')
  
  coordinates = np.array([source[n] for n in ['xt','yt','zt']]).T;
  #intensities = source['source'];
  
  #%% Unweighted 
  
  voxelization_parameter = dict(
        shape = io.shape(annotation_file), 
        dtype = None, 
        weights = None,
        method = 'sphere', 
        radius = (7,7,7), 
        kernel = None, 
        processes = None, 
        verbose = True
      )
  
  vox.voxelize(coordinates, sink=ws.filename('density', postfix='counts'), **voxelization_parameter);
  
  
  #%% 
  
  #p3d.plot(ws.filename('density', postfix='counts'))
  
  
  #%% Weighted 
  '''
  voxelization_parameter = dict(
        shape = io.shape(annotation_file),
        dtype = None, 
        weights = intensities,
        method = 'sphere', 
        radius = (7,7,7), 
        kernel = None, 
        processes = None, 
        verbose = True
      )
  
  vox.voxelize(coordinates, sink=ws.filename('density', postfix='intensities'), **voxelization_parameter);
  '''
  #%%
  
  #p3d.plot(ws.filename('density', postfix='intensities'))

  
#%% remove temp files 
#because apparently tempfile does not delete them by itself when used with Spyder
io.delete_file(ws.filename('stitched'))
io.delete_file(ws.filename('autofluo'))