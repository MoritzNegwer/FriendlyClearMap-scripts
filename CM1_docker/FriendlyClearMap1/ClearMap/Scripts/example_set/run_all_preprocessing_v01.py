#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 13:10:17 2020

@author: Moritz
"""
import os,numpy,sys


files = [
        '/home/wirrbel/ClearMap/ClearMap/Scripts/2020-09-02_SST_P25_Realigned+Arivis/S1_57645/process_template.py',
        '/home/wirrbel/ClearMap/ClearMap/Scripts/2020-09-02_SST_P25_Realigned+Arivis/S2_57646/process_template.py',
        '/home/wirrbel/ClearMap/ClearMap/Scripts/2020-09-02_SST_P25_Realigned+Arivis/S3_57647/process_template.py',
        '/home/wirrbel/ClearMap/ClearMap/Scripts/2020-09-02_SST_P25_Realigned+Arivis/S4_57648/process_template.py',
        '/home/wirrbel/ClearMap/ClearMap/Scripts/2020-09-02_SST_P25_Realigned+Arivis/S5_57654/process_template.py',
		'/home/wirrbel/ClearMap/ClearMap/Scripts/2020-09-02_SST_P25_Realigned+Arivis/S6_57657/process_template.py'
         ]


for brain in files:
    try:
        runfile(brain)
    except KeyError:
        pass
    
