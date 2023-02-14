# -*- coding: utf-8 -*-
"""
Template to run the processing pipeline
"""

#Load ClearMap Module
#TODO: Adjust path to your ClearMap folder 
runfile('/home/wirrbel/ClearMap/docs/conf.py') # Note that this is the only way I could get the ClearMap modules to load under Ubuntu 18.10
import ClearMap

import imagej
#load the parameters:
#TODO: adjust to point to your parameter file. 
execfile('/home/wirrbel/ClearMap/ClearMap/Scripts/2018-12-29_alignment_test/02_parameter_file.py')

