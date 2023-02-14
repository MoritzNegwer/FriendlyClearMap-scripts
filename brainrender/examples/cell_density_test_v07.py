"""
    This example shows how to use a PointsDensity
    actor to show the density of labelled cells
"""
# We begin by adding the current path to sys.path to make sure that the imports work correctly
import sys
sys.path.append('/home/wirrbel/brainrender_update_2022-01/')
import os
from shutil import copyfile
import in_place


import random
import numpy as np

from brainrender import Scene
from brainrender.actors import Points, PointsDensity

from rich import print
from myterial import orange
from pathlib import Path
import myterial



print(f"[{orange}]Running example: {Path(__file__).name}")

'''
def get_n_random_points_in_region(region, N):
    """
    Gets N random points inside (or on the surface) of a mes
    """

    region_bounds = region.mesh.bounds()
    X = np.random.randint(region_bounds[0], region_bounds[1], size=10000)
    Y = np.random.randint(region_bounds[2], region_bounds[3], size=10000)
    Z = np.random.randint(region_bounds[4], region_bounds[5], size=10000)
    pts = [[x, y, z] for x, y, z in zip(X, Y, Z)]

    ipts = region.mesh.insidePoints(pts).points()
    return np.vstack(random.choices(ipts, k=N))
'''

#hardcode paths 
# check for yourself by running "whereis transformix" in the command line 
transformix_bin = '/usr/bin/transformix'

def copy_and_optimize_transformParameters (new_folder,transform): 
    #takes the list of TransformParametersFile, copies them to the new folder and adapts their pointers 
    
    #copy over transform files 
    for transforms in transform:
        #
        file_to_copy = os.path.join(new_folder,os.path.split(transforms)[1])
        copyfile (transforms, file_to_copy)
        
        
        # adapt the InitialTransform parameter so that it points to the local file 
        # search and replace paths to InitialTransforms 
        with in_place.InPlace(file_to_copy) as file:
            for line in file:
                if ('Initial' in line) and ('NoInitialTransform' not in line):
                    #if 'NoInitialTransform' not in line:
                        #replace area between first quote and last slash with directory path
                        first_quote = line.find(' \"')
                        last_slash = line.rfind('/') 
                        line = line[:first_quote] + ' \"' + new_folder + line[last_slash:]
                file.write(line)
        file.close()
    

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

def transform_points (cellsfile, transform):
        #make temporary text file 
        cells_folder, file_name = os.path.split(cellsfile)
        #define filename with .txt extension
        txt_name = file_name[:-4] + ".txt"
        
        #make new subfolder 
        new_folder = os.path.join(cells_folder,'Aligned_CCF3')
        
        #if it doesn't exist yet, make new_folder 
        if not os.path.exists(new_folder):
            os.mkdir(new_folder)
        
        #save np file as txt 
        points = np.load(cellsfile)
        filename = os.path.join(new_folder,txt_name)
        
        with open(filename, 'w') as pointfile:
            #if indices:
            #    pointfile.write('index\n')
            #else:
            pointfile.write('point\n')
        
            pointfile.write(str(points.shape[0]) + '\n');
            np.savetxt(pointfile, points, delimiter = ' ', newline = '\n', fmt = '%.5e')
            pointfile.close();
            #np.savetxt(os.path.join(new_folder,txt_name), np.load(cellsfile), delimiter = ' ', newline = '\n', fmt = '%.5e')
        
        #copy transforms over and insert proper paths 
        copy_and_optimize_transformParameters (new_folder,transform)
        
        
        ##==== Step 1: intermediate alignment step to result.tif (translated atlas) === 
        
        #define copied transform 
        copied_transform = os.path.join(new_folder,os.path.split(transform[0])[1])
        
        #create transformix command
        intermediate_align = (transformix_bin
                     + ' -def ' + os.path.join(new_folder,txt_name) 
                     + ' -tp ' + copied_transform
                     + ' -out ' + new_folder)
        
        # run command 
        res = os.system(intermediate_align);
        
        ###=== Step 2: Re-save as new points file
        #load as np array 
        points_intermediate_align = parseElastixOutputPoints(os.path.join(new_folder,'outputpoints.txt'),indices=True)
        
        #write to intermediate points file 
        filename = os.path.join(new_folder,'transformed_points_intermediate.txt')
                
        with open(filename, 'w') as pointfile:
            pointfile.write('point\n')
            pointfile.write(str(points_intermediate_align.shape[0]) + '\n');
            np.savetxt(pointfile, points_intermediate_align, delimiter = ' ', newline = '\n', fmt = '%.5e')
            pointfile.close();
        
        ###=== Step 3: final (inverse) alignment 
        final_align = (transformix_bin
                       + ' -def ' + os.path.join(new_folder,'transformed_points_intermediate.txt')
                       + ' -tp ' +  os.path.join(new_folder,'TransformParameters.1.txt')
                       + ' -out ' + new_folder)
        
        res = os.system(final_align)
        
        ###=== Step 4: Load finally aligned cells 
        points_finally_aligned = parseElastixOutputPoints(os.path.join(new_folder,'outputpoints.txt'),indices=True)


        return points_finally_aligned 

def read_aligned_points(cellsfile):
    #subset of transform_points if you are really sure the transform has happened already
    cells_folder, file_name = os.path.split(cellsfile)
    
    #make new subfolder 
    new_folder = os.path.join(cells_folder,'Aligned_CCF3')
    
    points_finally_aligned = parseElastixOutputPoints(os.path.join(new_folder,'outputpoints.txt'),indices=True)

    return points_finally_aligned 



# Get a numpy array with (fake) coordinates of some labelled cells
#mos = scene.add_brain_region("MOs", alpha=0.0)
#coordinates = get_n_random_points_in_region(mos, 2000)

##DEBUG## 

target_mice = [["/media/wirrbel/Moritz_Negwer/SST_Maren/2019-04-23_SST_P25_Sample1_57645_Het/cells_transformed_to_Atlas_aligned.npy", "green"], #Het
               #["/media/wirrbel/Moritz_Negwer/SST_Maren/2019-04-24_SST_P25_Sample2_57646_WT/cells_transformed_to_Atlas_aligned.npy", "green"], #WT
               #["/media/wirrbel/Moritz_Negwer/SST_Maren/2019-04-25_SST_P25_Sample3_57647_WT/cells_transformed_to_Atlas_aligned.npy", "green"], #WT
               #["/media/wirrbel/Moritz_Negwer/SST_Maren/2019-04-25_SST_P25_Sample4_57648_Het/cells_transformed_to_Atlas_aligned.npy", "green"], #Het
               #["/media/wirrbel/Moritz_Negwer/SST_Maren/2019-04-26_SST_P25_Sample5_57654_WT/cells_transformed_to_Atlas_aligned.npy", "green"], #WT
               #["/media/wirrbel/Moritz_Negwer/SST_Maren/2019-04-26_SST_P25_Sample6_57657_Het/cells_transformed_to_Atlas_aligned.npy", "green"], #Het
               ]

target_folder = "/media/wirrbel/MN_2/2022-01-16_brainrender_tech_paper/SST_Maren/"

target_transformation = ["/home/wirrbel/2021-03-29_brainrender_2_preprocessing/elastix_files/Par_rotate90degaroundX_CCF3_nrrd_directions.txt",
                         "/home/wirrbel/2021-03-29_brainrender_2_preprocessing/Kim_ref_P21_v2_brain_CCF3/TransformParameters.0.txt", 
                         "/home/wirrbel/2021-03-29_brainrender_2_preprocessing/Kim_ref_P21_v2_brain_CCF3/TransformParameters.1.txt"]

target_regions = "Isocortex"

#### main function call 
for mouse in target_mice:
    #transform to CCF3 
    #cells = transform_points(mouse[0],target_transformation)
    #read transformed points (assuming all are transformed already)
    cells = read_aligned_points(mouse[0])
   # for region in target_regions:
        #render_screenshot (target_folder, cells, mouse[0], mouse[1], region, camera=cameras.get(region))
        #plane = render_screenshot (target_folder, cells, mouse[0], mouse[1], region, camera=techpaper_cam_01)
        
#multiply to get to µm value(assuming atlas has 25µm voxels)
cells = cells*25

# Add to scene
scene = Scene(screenshots_folder="/media/wirrbel/MN_2/2022-01-16_brainrender_tech_paper_densities/test")


##subset for region of interest 
#region to subset (added, but not displayed)
region_subset = scene.add_brain_region(target_regions, alpha=0.05,color=myterial.grey,hemisphere='left')
#subset for region 
cells = region_subset.mesh.insidePoints(cells).points()
 
scene.add(Points(cells, name="CELLS", colors="salmon",radius=2))

#slice one hemisphere 
plane = scene.atlas.get_plane(plane="sagittal", norm=[0, 0, -1])
#set sagittal cutting plane ever so slightly off the center 
#to prevent a render error with the "cartoon" render setting (see settings.py)
plane.mesh.SetPosition(7829.739797878185, 4296.026746612369, -5700.497969379508)
scene.slice(plane,close_actors=False)

#Changed points.py to be able to adapt color map, standard = "Dark2" from vedo: https://vedo.embl.es/autodocs/_modules/vedo/colors.html 
scene.add(PointsDensity(cells,colormap="twilight",radius=750))

# render
scene.render()
#scene.screenshot(name="test2")