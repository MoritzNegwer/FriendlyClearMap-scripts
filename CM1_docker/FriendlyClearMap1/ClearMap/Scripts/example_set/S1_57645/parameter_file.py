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
BaseDirectory = '/media/sf_E_DRIVE/SST_Maren/test_S1/';

#Data File and Reference channel File, usually as a sequence of files from the microscope
#Use \d{4} for 4 digits in the sequence for instance. As an example, if you have cfos-Z0001.ome.tif :
#os.path.join() is used to join the BaseDirectory path and the data paths:
# TODO Maren: 
cFosFile = os.path.join('/media/sf_E_DRIVE/SST_Maren/test_S1/190423_Sample1_57645_Het_TdTom_16-15-57/16-15-57_Sample1_57645_Het_TdTom_UltraII_C00_xyz-Table Z\d{4}.ome.tif');
AutofluoFile = cFosFile;

#convert csv to npy
cellspoints = numpy.loadtxt('/media/sf_E_DRIVE/SST_Maren/test_S1/cells.csv', delimiter = ",")
numpy.save(os.path.join(BaseDirectory,'cells.npy'), cellspoints)

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
PathReg        = '/media/sf_E_DRIVE/SST_Maren/atlas_p21/';
AtlasFile      = os.path.join(PathReg, 'Kim_ref_P21_v2_brain.tif');
LabeledImage   = os.path.join(PathReg, 'Kim_ref_P21_v2.1_label.nrrd'); 
AnnotationFile = os.path.join(PathReg, 'Regions_IDs_collapse.csv');


######################### Cell Detection Parameters using custom filters

#Spot detection method: faster, but optimised for spherical objects.
#You can also use "Ilastik" for more complex objects
ImageProcessingMethod = "SpotDetection";

#For illumination correction (necessitates a specific calibration curve for your microscope)
correctIlluminationParameter = {
    "flatfield"  : None,  # (True or None)  flat field intensities, if None do not correct image for illumination 
    "background" : None, # (None or array) background image as file name or array, if None background is assumed to be zero
    "scaling"    : None, # (str or None)        scale the corrected result by this factor, if 'max'/'mean' scale to keep max/mean invariant
    "save"       : None,       # (str or None)        save the corrected image to file
    "verbose"    : True    # (bool or int)        print / plot information about this step 
}

#Remove the background with morphological opening (optimised for spherical objects)
removeBackgroundParameter = {
    "size"    : (7,7),  # size in pixels (x,y) for the structure element of the morphological opening
    "save"    : os.path.join(BaseDirectory, 'background/background_7x7_\d{4}.tif'),     # file name to save result of this operation
    "save"    : None,
    "verbose" : True  # print / plot information about this step       
}

#Difference of Gaussians filter: to enhance the edges. Useful if the objects have a non smooth texture (eg: amyloid deposits)
filterDoGParameter = {
    "size"    : None,        # (tuple or None)      size for the DoG filter in pixels (x,y,z) if None, do not correct for any background
    "sigma"   : None,        # (tuple or None)      std of outer Gaussian, if None automatically determined from size
    "sigma2"  : None,        # (tuple or None)      std of inner Gaussian, if None automatically determined from size
    "save"    : None,        # (str or None)        file name to save result of this operation if None dont save to file 
    "verbose" : False      # (bool or int)        print / plot information about this step
}

#Extended maxima: if the peak intensity in the object is surrounded by smaller peaks: avoids overcounting "granular" looking objects
findExtendedMaximaParameter = {
    "hMax"      : None,            # (float or None)     h parameter (for instance 20) for the initial h-Max transform, if None, do not perform a h-max transform
    "size"      : 5,             # (tuple)             size for the structure element for the local maxima filter
    "threshold" : 0,        # (float or None)     include only maxima larger than a threshold, if None keep all local maxima
    "save"      : None,         # (str or None)       file name to save result of this operation if None dont save to file 
    "verbose"   : True       # (bool or int)       print / plot information about this step
}

#If no cell shape detection and the maximum intensity is not at the gravity center of the object, look for a peak intensity around the center of mass. 
findIntensityParameter = {
    "method" : 'Max',       # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
    "size"   : (4,4,4)      # (tuple)             size of the search box on which to perform the *method*
}

#Object volume detection. The object is painted by a watershed, until reaching the intensity threshold, based on the background subtracted image
detectCellShapeParameter = {
    "threshold" : None,     # (float or None)      threshold to determine mask. Pixels below this are background if None no mask is generated
    "save" : None,
    #"save"    : os.path.join(BaseDirectory, 'cellshape/cellshape\d{4}.tif'),        # (str or None)        file name to save result of this operation if None dont save to file 
    "verbose"   : True      # (bool or int)        print / plot information about this step if None take intensities at the given pixels
}


## Parameters for cell detection using spot detection algorithm 
detectSpotsParameter = {
    "correctIlluminationParameter" : correctIlluminationParameter,
    "removeBackgroundParameter"    : removeBackgroundParameter,
    "filterDoGParameter"           : filterDoGParameter,
    "findExtendedMaximaParameter"  : findExtendedMaximaParameter,
    "findIntensityParameter"       : findIntensityParameter,
    "detectCellShapeParameter"     : detectCellShapeParameter
}





#################### Heat map generation

##Voxelization: file name for the output:
VoxelizationFile = os.path.join(BaseDirectory, 'points_voxelized_aligned.tif');

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
    "processes": 4 
};


#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes. Be careful of the memory footprint of each process!
    "processes" : 2,
   
    #chunk sizes: number of planes processed at once
    #Default: 100,50,32
    "chunkSizeMax" : 40,
    "chunkSizeMin" : 20,
    "chunkOverlap" : 10,

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
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_cfos_to_auto'),
    "movingPoints" : None,
    "fixedPoints" : None,
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
SpotDetectionParameter = {
    "source" : cFosFile,
    "sink"   : (os.path.join(BaseDirectory, 'cells-allpoints_aligned.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints_aligned.npy')),
    "detectSpotsParameter" : detectSpotsParameter
};
SpotDetectionParameter = joinParameter(SpotDetectionParameter, cFosFileRange)

ImageProcessingParameter = joinParameter(StackProcessingParameter, SpotDetectionParameter);

FilteredCellsFile = (os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy'));

TransformedCellsFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas_aligned.npy')

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

