# -*- coding: utf-8 -*-
"""
Template to run the processing pipeline
"""

import os,numpy, sys

#pre-execute to make sure that Clearmap modules can be loaded

#generic solution
#basedir = os.path.abspath(os.path.join(os.path.dirname(__file__),'../../..'))

#hardcoded for docker
basedir = os.path.abspath('/CloudMap/FriendlyClearmap/clearmap/')

#insert into system folder
sys.path.insert(0, basedir)

import ClearMap

#load the parameters:

#original, with parameter file in same folder:
#parameter_filename = '/CloudMap/FriendlyClearmap/clearmap/Clearmap/Scripts/docker/parameter_file_Ilastik_template.py'

#with parameter file in the data folder
parameter_filename = '/CloudMap/Data/parameter_file_Ilastik_template.py'
exec(open(parameter_filename).read())


#resampling operations:
#######################
#resampling for the correction of stage movements during the acquisition between channels:
resampleData(**CorrectionResamplingParameterCfos);
resampleData(**CorrectionResamplingParameterAutoFluo);
#
#Downsampling for alignment to the Atlas:
resampleData(**RegistrationResamplingParameter);

