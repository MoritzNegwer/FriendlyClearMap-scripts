#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 22:07:08 2022

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

#hardcode paths 
# check for yourself by running "whereis transformix" in the command line 
transformix_bin = '/usr/bin/transformix'

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

def read_aligned_points(cellsfile):
    #subset of transform_points if you are really sure the transform has happened already
    cells_folder, file_name = os.path.split(cellsfile)
    
    #make new subfolder 
    new_folder = os.path.join(cells_folder,'Aligned_CCF3')
    
    points_finally_aligned = parseElastixOutputPoints(os.path.join(new_folder,'outputpoints.txt'),indices=True)

    return points_finally_aligned 

def render_screenshot (screenshots_folder, cells, cells_path ,cells_color,region_to_extract,camera=None):
    #try to make save folder 
    try:
        os.mkdir(screenshots_folder)
    except:
        pass
        
    # Create a scene
    scene = Scene (title=None,screenshots_folder=screenshots_folder,inset=None)
    #scene1.add_cells_from_file('/home/wirrbel/2020-03-02_Lynn_PV_P14/2020-03-02_Lynn_PV_P14_m2/cells_transformed_CCF.csv', color="red", radius=10,alpha=0.5)
        
    #for P56 alignment 
    #cellsfile = "/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m3/cells_transformed_to_Atlas_aligned.npy"
    #cells = np.load(cellsfile)
    cells = cells
    
    #define region name to extract 
    #region_to_extract = 'HIP'
    region_to_extract = region_to_extract
    
    #extract name of parent folder from cellsfile (assumed to be brain name)
    brain_name = str(region_to_extract) + "_" + os.path.basename(os.path.split(cells_path)[0]) 
        
    #multiply to get to µm value(assuming atlas has 25µm voxels)
    cells = cells*25
    #swap axes for yxz 
    #cells = cells[:,[1,0,2]]
  
    #for P56, ABA#2 
    #flip LR axis 
    cells[:,2] = cells[:,2] * -1
    cells[:,2] = cells[:,2] + 11400
      
    ##subset for region of interest 
    #region to subset (added, but not displayed)
    region_subset = scene.add_brain_region(region_to_extract, alpha=0.05,color=myterial.grey,hemisphere='left')
    #subset for region 
    cells = region_subset.mesh.insidePoints(cells).points()
    
    #add to scene 
    #scene.add(Points(cells,colors=cells_color,size=2))
    
    #slice one hemisphere 
    plane = scene.atlas.get_plane(plane="sagittal", norm=[0, 0, -1])
    #set sagittal cutting plane ever so slightly off the center 
    #to prevent a render error with the "cartoon" render setting (see settings.py)
    plane.mesh.SetPosition(7829.739797878185, 4296.026746612369, -5700.497969379508)
    scene.slice(plane,close_actors=False)
    
    #Changed points.py to be able to adapt color map, standard = "Dark2" from vedo: https://vedo.embl.es/autodocs/_modules/vedo/colors.html 
    scene.add(PointsDensity(cells,colormap="twilight",radius=750))
 
    #scene.render(interactive=True)

    
    #if region has a special camera angle attached to it, use this one. Othewise, use default (see brainrender settings)
    if camera is not None:
        scene.render(camera=camera, interactive=False,offscreen=True)
        scene.screenshot(name=brain_name)
        scene.close()
    else: 
        scene.screenshot(name=brain_name)
        scene.close()
    
    '''
    #debug
    scene.render()
    #periodically report camerra parameters
    import sched
    #scheduler = sched.scheduler(time.time,time.sleep)
    def periodically_report_cam_params (sc):
        try: 
            print(get_camera_params(scene=scene))
        except: 
            pass
        #scheduler.enter(5,1,periodically_report_cam_params, (sc))
        
    #scheduler.enter(5,1,periodically_report_cam_params, (sc))
    #scheduler.run()
    
    #debug purposes only (see below for test SST)
    return plane
    '''
###

#alternative for density
techpaper_cam_02 = {     
     "pos": (5538, -10429, -47712),
     "viewup": (0, -1, 0),
     "clippingRange": (33868, 52160),
     "focalPoint": (6888, 3571, -5717),
     "distance": 44288
     }

techpaper_cam_03 = {
     "pos": (-14632, -20614, -35939),
     "viewup": (0, -1, 0),
     "clippingRange": (27940, 61510),
     "focalPoint": (6888, 3571, -5717),
     "distance": 44288
     }

cameras = [
    techpaper_cam_02,
    techpaper_cam_03
    ]


target_regions = [#"grey", #All brain
                  "HIP", #Hippocampal region
                  "Isocortex", #Cortex proper
                  #"CNU", #Striatum + GP
                  "TH", #Thalamus +
                  #"CTXsp", #cortical subplate 
                  #"BLA", #BAsolateral Amygdala
                  #"MB", #Superior / inferior colliculus
                  "SSp-bfd", #S1 Barrel Field
                  "SS", #S1 
                  "AUD", #A1
                  "VIS" #V1
                  ]

#define empty cells array 
cells_all = np.empty((1,3))

##=== PV P56 WT === 
target_mice = [#["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m1/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               #["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m2/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               #["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m3/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               #["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m4/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               #["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m5/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               #["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m6/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               ["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m7/cells_transformed_to_Atlas_aligned.npy", "red"], #WT
               ["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m8/cells_transformed_to_Atlas_aligned.npy", "red"], #WT
               ["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-17_Lynn_PV_P50_m1/cells_transformed_to_Atlas_aligned.npy", "red"], #WT
               #["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-17_Lynn_PV_P50_m2/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               #["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m3/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               ["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m4/cells_transformed_to_Atlas_aligned.npy", "red"], #WT
               #["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m5/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               ["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m6/cells_transformed_to_Atlas_aligned.npy", "black"], #WT
               ["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-24_Lynn_PV_P50_m7/cells_transformed_to_Atlas_aligned.npy", "red"], #WT
               #["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-24_Lynn_PV_P50_m8/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               ]

target_folders = [
    #"/media/wirrbel/MN_2/2022-01-16_brainrender_tech_paper_densities/PV_P56_Averages_techpaper_cam_02/", #must end on trailing slash
    "/media/wirrbel/MN_2/2022-02-27_brainrender_for_kidsIQ/"
    ]

target_transformation = ["/home/wirrbel/2021-03-29_brainrender_2_preprocessing/elastix_files/Par_rotate90degaroundX_Clearmap1.txt"]

#### main function call 
for mouse in target_mice:
    #transform to CCF3 
    #cells = transform_points(mouse[0],target_transformation)
    #read transformed points (assuming all are transformed already)
    cells_temp = read_aligned_points(mouse[0])
    cells_all = np.concatenate((cells_temp,cells_all),axis=0)

#ugly hack for names 
name_path = "/media/wirrbel/MN_2/2022-01-16_brainrender_tech_paper_densities/names/PV_P56_WT_avg/"

for target_folder in target_folders:
    #once all added up, make 
    for region in target_regions:
        render_screenshot (target_folder, cells_all, name_path, myterial.black, region, camera=techpaper_cam_03)

#debug
#region = "Isocortex"
#render_screenshot (target_folder, cells_all, "PV_P56_WT_avg", myterial.red_dark, region, camera=techpaper_cam_02)

### === PV P56 Het === 

#define empty cells array 
cells_all = np.empty((1,3))

##PV P56 
target_mice = [["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m1/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               ["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m2/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               ["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m3/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               ["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m4/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               ["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m5/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               ["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m6/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               #["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m7/cells_transformed_to_Atlas_aligned.npy", "red"], #WT
               #["/media/wirrbel/MN_2/Lynn_PV_P50/Lukas_PV_50_m8/cells_transformed_to_Atlas_aligned.npy", "red"], #WT
               #["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-17_Lynn_PV_P50_m1/cells_transformed_to_Atlas_aligned.npy", "red"], #WT
               ["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-17_Lynn_PV_P50_m2/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               ["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m3/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               #["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m4/cells_transformed_to_Atlas_aligned.npy", "red"], #WT
               ["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m5/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               #["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-19_Lynn_PV_P50_m6/cells_transformed_to_Atlas_aligned.npy", "black"], #WT
               #["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-24_Lynn_PV_P50_m7/cells_transformed_to_Atlas_aligned.npy", "red"], #WT
               ["/media/wirrbel/MN_2/Lynn_PV_P50/2020-02-24_Lynn_PV_P50_m8/cells_transformed_to_Atlas_aligned.npy", "red"], #Het
               ]

target_folders = [
    #"/media/wirrbel/MN_2/2022-01-16_brainrender_tech_paper_densities/PV_P56_Averages_techpaper_cam_02/", #must end on trailing slash
    "/media/wirrbel/MN_2/2022-02-27_brainrender_for_kidsIQ/"
    ]

target_transformation = ["/home/wirrbel/2021-03-29_brainrender_2_preprocessing/elastix_files/Par_rotate90degaroundX_Clearmap1.txt"]

#### main function call 
for mouse in target_mice:
    #transform to CCF3 
    #cells = transform_points(mouse[0],target_transformation)
    #read transformed points (assuming all are transformed already)
    cells_temp = read_aligned_points(mouse[0])
    cells_all = np.concatenate((cells_temp,cells_all),axis=0)

#ugly hack for names 
name_path = "/media/wirrbel/MN_2/2022-01-16_brainrender_tech_paper_densities/names/PV_P56_Het_avg/"

for target_folder in target_folders:
    #once all added up, make 
    for region in target_regions:
        render_screenshot (target_folder, cells_all, name_path, myterial.red_dark, region, camera=techpaper_cam_03)

#debug
#region = "Isocortex"
#render_screenshot (target_folder, cells_all, "PV_P56_WT_avg", myterial.red_dark, region, camera=techpaper_cam_02)
