# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 16:17:06 2020

@author: AnniePutman
"""

#write some numpy serialized data files 
import numpy as np
import os
import pandas as pd

dataloc = r"E:\WorkComputer\Documents\Documents\MultidecadalChangesPackage\Data"
os.chdir(dataloc)
variables = {'Lons':Lons, 'Lats':Lats, 'icyears':icyears, 'd18O':d18O, 'Pre':Pre, 
             'Tmp':Tmp, 'PW':PW, 'grids':grids, 
     'pregrids':pregrids, 'tmpgrids': tmpgrids, 'pwgrids':pwgrids, 'glats':glats, 
     'glons':glons, 'cgseas':cgseas, 'cgpre':cgpre, 'pretot':pretot, 'cgnmons':cgnmons, 
     'cgtcov':cgtcov, 'years':years, 'mons':mons, 'cgtmp':cgtmp, 'cgprjs':cgprjs}



for variable in variables.keys():
    if not os.path.exists(variable):
        if variables[variable] is not np.array:
            variables[variable] = np.array(variables[variable])
        np.save(variable, variables[variable])

df.to_csv('df.csv')
