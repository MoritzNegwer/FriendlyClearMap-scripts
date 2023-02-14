#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 02:13:31 2022

@author: wirrbel
"""

# We begin by adding the current path to sys.path to make sure that the imports work correctly
import sys
sys.path.append('/home/wirrbel/brainrender_update_2022-01/')
import os
import pandas as pd
import numpy as np
from shutil import copyfile
import in_place

from vedo import embedWindow  # for more explanations about these two lines checks the notebooks workflow example
embedWindow(None)

# Import variables
from brainrender import * # <- these can be changed to personalize the look of your renders

# Import brainrender classes and useful functions
from brainrender import Scene
from brainrender.actors import Points, PointsDensity, Volume

#for colors
import myterial

import tifffile
import matplotlib.pyplot as plt 


def parseElastixOutputPoints(filename, indices = True):
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
        return numpy.zeros((0,3));

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

def simple_density(cells_df, target_dir, stack_template,mouse_id):

    #create an empty stack with the same shape as the template 
    stack_dens = np.zeros(shape=stack_template.shape).astype('uint16')
    
    #calculate density per voxel (only for cortex)
    density_per_voxel = cells_df.loc[cells_df['is_cortex']][['x','y','z']].value_counts()
    
    #put into empty stack 
    for voxel in density_per_voxel.iteritems():
        coords = voxel[0]
        try:
            stack_dens[coords[2],coords[0],coords[1]] = voxel[1]
        except IndexError:
            #if the cell is outside of the boundaries, disregard
            pass
        
    #write to image stack 
    tifffile.imwrite(os.path.join(target_dir,mouse_id+'_density_stack_v01.tif'),stack_dens)

def up_one(area_id,region_ids):
    return int(region_ids.loc[area_id]['parent_id'])

def cortex_or_not(area_id,region_ids,target_id=500,target_name='Isocortex'):
    while region_ids.loc[area_id]['structure_order'] > target_id:
        #still lower than cortex, move one step up 
        area_id = up_one(area_id,region_ids)
        
    #once the loop is broken, the chain of parent ids has either led to the cortex or elsewhere
    if region_ids.loc[area_id]['name'] == target_name:
        return True
    else:
        return False

def cortical_depth_map(cells_dir,cortical_distances_file,annotations_file,regions_ids,genotype,age,target_folder):
    #load CCF3-transformed cells
    cells = parseElastixOutputPoints(os.path.join(cells_dir,'Aligned_CCF3/outputpoints.txt')).astype('uint16')
    
    #load cortical distance (from ClearMap2 resources)
    distances = tifffile.imread(cortical_distances_file).astype('float32')
    
    #determine distance to surface (i.e. cortical depth) for each cell 
    cell_values = []
    for cell in cells:
        try:
            value = distances[cell[2],cell[0],cell[1]]
            cell_values.append(value)
        except IndexError:
            #if there is an error, the cell coordinates are outside the atlas. Give distance as 0 (needed for consistency reasons, cell will be discarded anyways)
            cell_values.append(0)
            pass
        
    #read atlas and assign region to each 
    ccf3 = tifffile.imread(annotations_file).astype('uint32')
    cell_areas = []
    for cell in cells:
        try:
            area = ccf3[cell[2],cell[0],cell[1]]
            cell_areas.append(area)
        except IndexError:
            #if there is an error, the cell coordinates are outside the atlas. Count as background. 
            cell_areas.append(0)
        
    #load region id
    region_ids = pd.read_csv(regions_ids,index_col='id')
    
    #pull everything together in a pandas dataframe 
    cells_df = pd.DataFrame(data=cells,columns={'x','y','z'})
    cells_df = cells_df.assign(structure_id='',structure_name='',structure_acronym='',structure_order='',parent_id='',mouse='',genotype='',age='')
    cells_df = cells_df.assign(structure_id = cell_areas)
    cells_df = cells_df.assign(is_cortex='')
    cells_df = cells_df.assign(pia_distance=cell_values)
        
    #generate mouse ID 
    mouse_id = os.path.basename(os.path.dirname(cells_dir))
    cells_df['mouse'] = mouse_id
    cells_df['genotype'] = genotype
    cells_df['age'] = age
    
    
    #check for every cell whether it sits in cortex or not, add in dataframe 
    for index,cell in cells_df.iterrows():
        try: 
            area_id = cell['structure_id']
            is_ctx = cortex_or_not(area_id,region_ids)
            #when encountering areas not in the list (happens with some obscurer ones)
        except KeyError:
            #print('skipping area: ', area_id)
            is_ctx = False
        
        if is_ctx:
            cells_df.at[index,'is_cortex'] = True
        else:
            cells_df.at[index,'is_cortex'] = False
    
    #create a cortex-subset 
    cortex_subset = cells_df.loc[cells_df['is_cortex']]
    
    #define square fig 
    plt.figure(figsize = (10,10), dpi = 300)
    #plot a histogram
    plt.hist(cortex_subset['pia_distance'].to_numpy()*25,bins=40,range={0,1600},histtype='step',orientation='horizontal')
    #invert y-axis
    plt.gca().invert_yaxis()
    #set maximum x limit
    plt.xlim(0,20000)
    #label x+y axes
    plt.xlabel('Cells')
    plt.ylabel('Cortical depth (µm)')
    #save the histogram
    plt.savefig(os.path.join(target_folder,mouse_id+'_cortical_depth_v01.svg'))
    
    #save cells_df as pickle
    cells_df.to_pickle(os.path.join(target_folder,mouse_id+'_v01.pickle'))
    #save cells_df as csv
    cells_df.to_csv(os.path.join(target_folder,mouse_id+'_v01.csv'))
    
    #save heatmap tiff stack 
    simple_density(cells_df,target_folder,ccf3,mouse_id)
    
    return mouse_id, cortex_subset

def draw_all_depth_maps(all_depths_file,output_folder):
    #standalone file for creating histograms of the collected files 
    #read pickle file 
    all_depths = pd.read_pickle(all_depths_file)
    
    #when all_depths is missing, manually assemble from individual pickles
    '''
    depth_pickles = glob.glob('/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/2022-10-13_cortical_depth_v01/*.pickle')
    all_mice_depth = pd.DataFrame()
    for pickle_file in depth_pickles:
        #split filename off
        pathname, filename = os.path.split(pickle_file)
        #extract mouse ID from filename
        mouse_id = os.path.splitext(filename)[0]
        #load the pickle into a Pandas dataframe
        temp_df = pd.read_pickle(pickle_file)
        #extract the cortical depth markers from temp_df
        cortical_depth_temp = temp_df.loc[temp_df['is_cortex']]['pia_distance']
        cortical_depth_temp = cortical_depth_temp.rename(mouse_id)
        #put into all_mice_depth
        all_mice_depth = pd.concat([all_mice_depth, cortical_depth_temp], axis=1)
    '''
    
    #genotype_list
    PV_P14_WT = ['2020-03-02_Lynn_PV_P14_m1',
     '2020-03-02_Lynn_PV_P14_m2',
     '2020-03-02_Lynn_PV_P14_m3',
     '2020-03-02_Lynn_PV_P14_m5',
     'Lukas_PV_P14_m1_66632',
     'Lukas_PV_P14_m2_66633',
     'Lukas_PV_P14_m6_190926']
    
    PV_P14_Het = ['2020-03-02_Lynn_PV_P14_m4',
     'Lukas_PV_P14_m3_66634',
     'Lukas_PV_P14_m4_66635',
     'Lukas_PV_P14_m5_190925']
    
    PV_P50_WT = ['Lukas_PV_50_m6',
     'Lukas_PV_50_m7',
     'Lukas_PV_50_m8',
     '2020-02-17_Lynn_PV_P50_m1',
     '2020-02-19_Lynn_PV_P50_m4',
     '2020-02-19_Lynn_PV_P50_m6',
     '2020-02-24_Lynn_PV_P50_m7']
    
    PV_P50_Het = ['Lukas_PV_50_m1',
     'Lukas_PV_50_m2',
     'Lukas_PV_50_m3',
     'Lukas_PV_50_m4',
     'Lukas_PV_50_m5',
     '2020-02-17_Lynn_PV_P50_m2',
     '2020-02-19_Lynn_PV_P50_m3',
     '2020-02-19_Lynn_PV_P50_m5',
     '2020-02-24_Lynn_PV_P50_m8']
    
    SST_WT = ['2019-04-24_SST_P25_Sample2_57646_WT',
    '2019-04-25_SST_P25_Sample3_57647_WT',
    '2019-04-26_SST_P25_Sample5_57654_WT']
    
    SST_Het = ['2019-04-23_SST_P25_Sample1_57645_Het',
    '2019-04-25_SST_P25_Sample4_57648_Het',
    '2019-04-26_SST_P25_Sample6_57657_Het']
    
    VIP_WT = ['Lynn_VIP_P26_m69366_WT',
      'VIP_P26_m69369_WT_tdtom',
     'Lynn_VIP_P26_m69370_WT',
     'Lynn_VIP_P26_m69371_WT',
     'm1',
     'm2',
     'm3']
    
    VIP_Het = ['Lynn_VIP_P26_m69368_Het',
     'm4',
     'm5',
     'm6']

    genotypes_list = [PV_P14_WT,PV_P14_Het,PV_P50_WT,PV_P50_Het,SST_WT,SST_Het,VIP_WT,VIP_Het]

    #create one plot for each line 
    for genotype in genotypes_list:
        #when manually assembled from individual pickles, append _v01 to every element in genotypes_list
        #genotype = [s + '_v01' for s in genotype]
            
        #get data if run in-line
        subset = all_depths[genotype]
        #get first mouseID as proxy for genotype
        genotype_id = subset.columns.to_list()[0]
        weights = None
        
        #optional: If you want to make a single line for all mice
        number_mice = len(subset.columns.to_list())
        subset = subset.stack()
        #then generate weights to normalize for the number of mice 
        weights = np.ones(subset.shape[0])/number_mice

        #make figure
        #define square fig 
        plt.figure(figsize = (10,10), dpi = 300)
        #plot a histogram
        plt.hist(subset.to_numpy()*25,bins=40,weights=weights,range={0,1600},histtype='step',orientation='horizontal')
        #invert y-axis
        plt.gca().invert_yaxis()
        #set maximum x limit
        plt.xlim(0,15000)
        #label x+y axes
        plt.xlabel('Cells')
        plt.ylabel('Cortical depth (µm)')
        #save the histogram
        plt.savefig(os.path.join(target_folder,'all_' + genotype_id +'_cortical_depth_v01.svg'))
        

cortical_distances_file = '/home/wirrbel/ClearMap2/ClearMap/Resources/Atlas/ABA_25um_distance_to_surface.tif'
annotations_file = '/home/wirrbel/ClearMap2/ClearMap/Resources/Atlas/ABA_25um_annotation.tif'
regions_ids = '/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_docker/atlas_p21/regions_IDs.csv'
target_folder = '/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/2022-10-13_cortical_depth_v01/'

### from tsne_plot -> assign_cells_to_regions_v01

all_depths = pd.DataFrame()
#==== P14 ====
WorkFoldersP14 = [['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P14/2020-03-02_Lynn_PV_P14_m1/','WT',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P14/2020-03-02_Lynn_PV_P14_m2/','WT',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P14/2020-03-02_Lynn_PV_P14_m3/','WT',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P14/2020-03-02_Lynn_PV_P14_m4/','Het',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P14/2020-03-02_Lynn_PV_P14_m5/','WT',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lukas_PV_P14/Lukas_PV_P14_m1_66632/','WT',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lukas_PV_P14/Lukas_PV_P14_m2_66633/','WT',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lukas_PV_P14/Lukas_PV_P14_m3_66634/','Het',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lukas_PV_P14/Lukas_PV_P14_m4_66635/','Het',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lukas_PV_P14/Lukas_PV_P14_m5_190925/','Het',14],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lukas_PV_P14/Lukas_PV_P14_m6_190926/','WT',14]]


for brain in WorkFoldersP14:
   
    #disassemble into components 
    cells_dir = brain[0]
    genotype = brain[1]
    age = brain[2]
    
    mouse_id, cortex_subset = cortical_depth_map(cells_dir,cortical_distances_file,annotations_file,regions_ids,genotype,age,target_folder)
    
    #add to all_depths
    all_depths = pd.concat([all_depths, cortex_subset['pia_distance'].rename(mouse_id)],axis=1)

    
#==== P56 ==== 
WorkFoldersP56 = [['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/Lukas_PV_50_m1/','Het',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/Lukas_PV_50_m2/','Het',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/Lukas_PV_50_m3/','Het',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/Lukas_PV_50_m4/','Het',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/Lukas_PV_50_m5/','Het',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/Lukas_PV_50_m6/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/Lukas_PV_50_m7/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/Lukas_PV_50_m8/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/2020-02-17_Lynn_PV_P50_m1/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/2020-02-17_Lynn_PV_P50_m2/','Het',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m3/','Het',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m4/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m5/','Het',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m6/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/2020-02-24_Lynn_PV_P50_m7/','WT',56],
               ['/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/Lynn_PV_P50/2020-02-24_Lynn_PV_P50_m8/','Het',56]
               ]

for brain in WorkFoldersP56:
   
    #disassemble into components 
    cells_dir = brain[0]
    genotype = brain[1]
    age = brain[2]
    
    mouse_id, cortex_subset = cortical_depth_map(cells_dir,cortical_distances_file,annotations_file,regions_ids,genotype,age,target_folder)
        
    #add to all_depths
    all_depths = pd.concat([all_depths, cortex_subset['pia_distance'].rename(mouse_id)],axis=1)
    
#===SST===
WorkFoldersSST = [["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/SST_Maren/2019-04-23_SST_P25_Sample1_57645_Het/", 'Het',28], #Het
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/SST_Maren/2019-04-24_SST_P25_Sample2_57646_WT/", 'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/SST_Maren/2019-04-25_SST_P25_Sample3_57647_WT/", 'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/SST_Maren/2019-04-25_SST_P25_Sample4_57648_Het/", 'Het',28], #Het
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/SST_Maren/2019-04-26_SST_P25_Sample5_57654_WT/", 'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/SST_Maren/2019-04-26_SST_P25_Sample6_57657_Het/", 'Het',28], #Het
               ]

for brain in WorkFoldersSST:
   
    #disassemble into components 
    cells_dir = brain[0]
    genotype = brain[1]
    age = brain[2]
    
    mouse_id, cortex_subset = cortical_depth_map(cells_dir,cortical_distances_file,annotations_file,regions_ids,genotype,age,target_folder)
    
    #add to all_depths
    all_depths = pd.concat([all_depths, cortex_subset['pia_distance'].rename(mouse_id)],axis=1)    
    
#===VIP===
WorkFoldersVIP = [["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_docker/Lynn_VIP_P26_m69366_WT/", 'WT', 28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_docker/Lynn_VIP_P26_m69368_Het/", 'Het', 28], #Het
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_docker/Lynn_VIP_P26_m69369_WT/VIP_P26_m69369_WT_tdtom/",'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_docker/Lynn_VIP_P26_m69370_WT/",'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_docker/Lynn_VIP_P26_m69371_WT/", 'WT', 28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_Arivis/m1/", 'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_Arivis/m2/", 'Het', 28], #Het
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_Arivis/m3/", 'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_Arivis/m4/", 'Het',28], #Het
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_Arivis/m5/", 'WT',28], #WT
               ["/media/wirrbel/MN_4/Friendly_Clearmap_Docker_All_Data_For_Tech_Paper/VIP_Arivis/m6/", 'Het',28], #Het
                ]

for brain in WorkFoldersVIP:
   
    #disassemble into components 
    cells_dir = brain[0]
    genotype = brain[1]
    age = brain[2]
    
    mouse_id, cortex_subset = cortical_depth_map(cells_dir,cortical_distances_file,annotations_file,regions_ids,genotype,age,target_folder)
    
    #add to all_depths
    all_depths = pd.concat([all_depths, cortex_subset['pia_distance'].rename(mouse_id)],axis=1)    


#=== END === 
#save cells_df as pickle
all_depths.to_pickle(os.path.join(target_folder,'all_depths_v01.pickle'))

#save cells_df as csv
all_depths.to_csv(os.path.join(target_folder,'all_depths_v01.csv'))
