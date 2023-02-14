# -*- coding: utf-8 -*-
"""
Example script to set up the parameters for the image processing pipeline
"""

######################### Import modules

import os, numpy, math

import ClearMap.Settings as settings
import ClearMap.IO as io

from ClearMap.Alignment.Resampling import resampleData;
from ClearMap.Alignment.Elastix import alignData, transformPoints
from ClearMap.ImageProcessing.CellDetection import detectCells
from ClearMap.Alignment.Resampling import resamplePoints, resamplePointsInverse
from ClearMap.Analysis.Label import countPointsInRegions
from ClearMap.Analysis.Voxelization import voxelize
from ClearMap.Analysis.Statistics import thresholdPoints
from ClearMap.Utils.ParameterTools import joinParameter
from ClearMap.Analysis.Label import labelToName


######################### Data parameters

#Directory to save all the results, usually containing the data for one sample
#Directory to save all the results, usually containing the data for one sample
BaseDirectory = '/CloudMap/Data/';

#Data File and Reference channel File, usually as a sequence of files from the microscope
#Use \d{4} for 4 digits in the sequence for instance. As an example, if you have cfos-Z0001.ome.tif :
#os.path.join() is used to join the BaseDirectory path and the data paths:
# TODO Maren: 

#AutofluoFile = (os.path.join(BaseDirectory,'Source/autofluo/15-45-35_2020-11-01_WFA_Mo_S01_autofluo_Horizontal_focus_UltraII_C00_xyz-Table Z\d{4}.ome.tif'))

cFosFile = (os.path.join(BaseDirectory,'Source/protein/17-11-30_m359920_tdtom_UltraII_C00_xyz-Table Z\d{4}.ome.tif');
AutofluoFile = cFosFile

#Specify the range for the cell detection. This doesn't affect the resampling and registration operations
cFosFileRange = {'x' : all, 'y' : all, 'z' : all};

#Resolution of the Raw Data (in um / pixel)
OriginalResolution = (3.25, 3.25, 3);

#Orientation: 1,2,3 means the same orientation as the reference and atlas files.
#Flip axis with - sign (eg. (-1,2,3) flips x). 3D Rotate by swapping numbers. (eg. (2,1,3) swaps x and y)
FinalOrientation = (1,2,3);

#Resolution of the Atlas (in um/ pixel)
AtlasResolution = (25, 25, 25);

#Path to registration parameters and atlases
PathReg        = '/Cloudmap/FriendlyClearmap/atlas/';
AtlasFile      = os.path.join(PathReg, 'template_25.tif');
LabeledImage   = os.path.join(PathReg, 'annotation_25_full.nrrd'); 
AnnotationFile = os.path.join(PathReg, 'regions_IDs.csv');



######################### Cell Detection Parameters using custom filters

#Spot detection method: faster, but optimised for spherical objects.
#You can also use "Ilastik" for more complex objects
ImageProcessingMethod = "Ilastik";

ilastikParameter = {
    "classifier" : "/CloudMap/classifiers/classifier_current.ilp",
    "classindex" : 0,
    "save"       : None, # (str or None)       file name to save result of this operation if None dont save to file 
    "verbose"    : True       # (bool or int)       print / plot information about this step
}


#If the maximum instensity is not at the gravity center of the object, look for a peak intensity around the center of mass. 
findCellIntensityParameter = {
    "method" : 'Max',       # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
    "verbose": True         # (bool or int)       print / plot information about this step
}


## Parameters for cell detection using spot detection algorithm 
detectCellsParameter = {
    "classifyCellsParameter"  : ilastikParameter,
    "findCellIntensityParameter"  : findCellIntensityParameter
}





#################### Heat map generation

##Voxelization: file name for the output:
VoxelizationFile = os.path.join(BaseDirectory, 'points_voxelized.tif');

# Parameter to calculate the density of the voxelization
voxelizeParameter = {
    #Method to voxelize
    "method" : 'Spherical', # Spherical,'Rectangular, Gaussian'
       
    # Define bounds of the volume to be voxelized in pixels
    "size" : (15,15,15),  

    # Voxelization weigths (e/g intensities)
    "weights" : None
};










############################ Config parameters

#Processes to use for Resampling (usually twice the number of physical processors)
ResamplingParameter = {
    "processes": 12 
};


#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes. Be careful of the memory footprint of each process!
    "processes" : 12,
   
    #chunk sizes: number of planes processed at once
    #optimized for 64GB RAM, can be higher if more available
    "chunkSizeMax" : 400,
    "chunkSizeMin" : 100,
    "chunkOverlap" : 20,

    #optimize chunk size and number to number of processes to limit the number of cycles
    "chunkOptimization" : True,
    
    #increase chunk size for optimization (True, False or all = automatic)
    "chunkOptimizationSize" : all,
   
    "processMethod" : "parallel"
   };






######################## Run Parameters, usually you don't need to change those


### Resample Fluorescent and CFos images
# Autofluorescent cFos resampling for aquisition correction

ResolutionAffineCFosAutoFluo =  (16, 16, 16);

CorrectionResamplingParameterCfos = ResamplingParameter.copy();

CorrectionResamplingParameterCfos["source"] = cFosFile;
CorrectionResamplingParameterCfos["sink"]   = os.path.join(BaseDirectory, 'cfos_resampled.tif');
    
CorrectionResamplingParameterCfos["resolutionSource"] = OriginalResolution;
CorrectionResamplingParameterCfos["resolutionSink"]   = ResolutionAffineCFosAutoFluo;

CorrectionResamplingParameterCfos["orientation"] = FinalOrientation;
   
   
   
#Files for Auto-fluorescence for acquisition movements correction
CorrectionResamplingParameterAutoFluo = CorrectionResamplingParameterCfos.copy();
CorrectionResamplingParameterAutoFluo["source"] = AutofluoFile;
CorrectionResamplingParameterAutoFluo["sink"]   = os.path.join(BaseDirectory, 'autofluo_for_cfos_resampled.tif');
   
#Files for Auto-fluorescence (Atlas Registration)
RegistrationResamplingParameter = CorrectionResamplingParameterAutoFluo.copy();
RegistrationResamplingParameter["sink"]            =  os.path.join(BaseDirectory, 'autofluo_resampled.tif');
RegistrationResamplingParameter["resolutionSink"]  = AtlasResolution;
   

### Align cFos and Autofluo

CorrectionAlignmentParameter = {            
    #moving and reference images
    "movingImage" : os.path.join(BaseDirectory, 'autofluo_for_cfos_resampled.tif'),
    "fixedImage"  : os.path.join(BaseDirectory, 'cfos_resampled.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(PathReg, 'Par0000affine_acquisition.txt'),
    "bSplineParameterFile" : None,
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_cfos_to_auto')
    }; 
  

### Align Autofluo and Atlas

#directory of the alignment result
RegistrationAlignmentParameter = CorrectionAlignmentParameter.copy();

RegistrationAlignmentParameter["resultDirectory"] = os.path.join(BaseDirectory, 'elastix_auto_to_atlas');
    
#moving and reference images
RegistrationAlignmentParameter["movingImage"]  = AtlasFile;
RegistrationAlignmentParameter["fixedImage"]   = os.path.join(BaseDirectory, 'autofluo_resampled.tif');

#elastix parameter files for alignment
RegistrationAlignmentParameter["affineParameterFile"]  = os.path.join(PathReg, 'Par0000affine.txt');
RegistrationAlignmentParameter["bSplineParameterFile"] = os.path.join(PathReg, 'Par0000bspline.txt');
RegistrationAlignmentParameter["movingPoints"] = os.path.join(BaseDirectory, 'atlas_landmarks.txt');
RegistrationAlignmentParameter["fixedPoints"] = os.path.join(BaseDirectory, 'autofluo_landmarks.txt');



# result files for cell coordinates (csv, vtk or ims)
CellDetectionParameter = {
    "method" : ImageProcessingMethod,
    "source" : cFosFile,
    "sink"   : (os.path.join(BaseDirectory, 'cells-allpoints_aligned.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints_aligned.npy')),
};
CellDetectionParameter = joinParameter(CellDetectionParameter, detectCellsParameter)
CellDetectionParameter = joinParameter(CellDetectionParameter, cFosFileRange)

ImageProcessingParameter = joinParameter(StackProcessingParameter, CellDetectionParameter);

FilteredCellsFile = (os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy'));

TransformedCellsFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas_alined.npy')

### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to referenece image size

CorrectionResamplingPointsParameter = CorrectionResamplingParameterCfos.copy();
CorrectionResamplingPointsParameter["pointSource"] = os.path.join(BaseDirectory, 'cells.npy');
CorrectionResamplingPointsParameter["dataSizeSource"] = cFosFile;
CorrectionResamplingPointsParameter["pointSink"]  = None;

CorrectionResamplingPointsInverseParameter = CorrectionResamplingPointsParameter.copy();
CorrectionResamplingPointsInverseParameter["dataSizeSource"] = cFosFile;
CorrectionResamplingPointsInverseParameter["pointSink"]  = None;

## Transform points from corrected to registered
# downscale points to referenece image size
RegistrationResamplingPointParameter = RegistrationResamplingParameter.copy();
RegistrationResamplingPointParameter["dataSizeSource"] = cFosFile;
RegistrationResamplingPointParameter["pointSink"]  = None;

