# -*- coding: utf-8 -*-
"""
Script that does all processing, analysis, makes figures, and produces statistics 
for the paper. 
"""
import os
import numpy as np
import importlib

rootdir = """E:\WorkComputer\Documents\Documents\Python\WaterIsotopeAnalysis_012016\GeospatialInterpolation\Manuscript_Clean"""
os.chdir(os.path.join(rootdir, 'Scripts'))
#this is the set of functions for getting iCAM and measurement data
import GetMSData as GetData
#this is the set of functions for making various figures
import MakeMSFigs as MakeFigures 
#get the slope and difference for data, model subset and full model
import GetMSChange as GetChange
#estimate residence time for iCAM full
import EstimateMSRT as EstimateRT 

import ProcessMSData as ProcessData

import GetMSStats as GetStats

#LOCATIONS OF THINGS

aggplots = os.path.join(rootdir, 'Figs')
ncfiles = os.path.join(rootdir, 'Data')
teleloc = os.path.join(rootdir, 'Data')
prefileindir = os.path.join(rootdir, 'Data\Reanalysis')
CAMdir = r"""E:\WorkComputer\WDrive\iCAM\forAnnie"""  #this is big so needs to be held elsewhere (now it can potentially be in Data)
getregrid = False
#%% ######## GET AND PROCESS DATA (MODEL AND MEASUREMENT) ##################

#Get and process model data (3 ensemble members)
(Lons, Lats, icyears, d18O, Pre, Tmp, PW, grids, 
     pregrids, tmpgrids, pwgrids) = ProcessData.processCAM(GetData, CAMdir)

#Get and process measurement data
(glats, glons, cgseas, cgpre, pretot, cgnmons, 
     cgtcov, years, mons, cgtmp, cgprjs) = ProcessData.processData(GetData, ncfiles)

startind = np.where(np.unique(icyears) == np.min(years))[0][0]
endind = np.where(np.unique(years) == np.max(icyears))[0][0]+1
gridsub = np.where(np.isnan(cgseas[:, :, :endind, :]), cgseas[:, :, :endind, :], grids[:, :, startind:, :])
MODSUBSurfEst, MODSUBSurfEstSE = EstimateRT.getSurfEst(gridsub)
#use the same time frame
MODSurfEst, MODSurfEstSE = EstimateRT.getSurfEst(grids[:, :, startind:, :])
#%% ######### RESIDENCE TIME ESTIMATE OF D18O TIMESERIES (MODEL AND MEASUREMENT) ################

#estimate residence time d18O timeseries from model data
d18O_ts, Dd18O_ts, yearrange = EstimateRT.EstModRT(icyears, grids, pregrids, pwgrids)

#estimate residence time for measurement from daily scale reanalysis data
SurfEst, SurfEstSE = EstimateRT.getSurfEst(cgseas)
if getregrid:
    EstimateRT.GetRTData(prefileindir)
    EstimateRT.regridRTdata()
#this is returning a lot of nans... something is wrong
Tseas, d18O_est, Dd18O_est = EstimateRT.getMeasRT(SurfEst, SurfEstSE, prefileindir, ncfiles)
cgseasRT = np.swapaxes(np.swapaxes(Dd18O_est, 2, 0), 3, 1)

#%% ############# DO REGRESSIONS OF MODEL, DATA, AND RESIDENCE TIME TIMESERIES ########
#These should come out in per mil per decade
#data
regresults = GetChange.GetMultidecadalChange(glats, glons, years, cgseas,                                              
                                             cgtcov, cgnmons)
#model d18O
MODregresults = GetChange.getCAMChange(icyears, years, Lats, 
                             Lons, grids, cgnmons, cgtcov, [0, 0])
#subset
MODSUBregresults = GetChange.GetMultidecadalChange(glats, 
                                             glons, years[years <= np.max(icyears)],  
                                             gridsub, cgtcov, cgnmons)

#model residence time
RTMODregresults = GetChange.getCAMChange(yearrange, years, Lats, Lons, 
                                         Dd18O_ts, cgtcov, cgnmons, [0,0])

#reanalysis residence time
RTregresults = GetChange.GetMultidecadalChange(glats, glons, years, cgseasRT,                                              
                                             cgtcov, cgnmons)
#summarize comparison of residence time estimate and model or measurement slope estimate
latareas, totalarea = GetStats.areaPerGridCell(glats, glons)
GetStats.summary(-RTregresults, regresults, latareas, totalarea)
GetStats.summary(-RTMODregresults, MODregresults, latareas, totalarea)
GetStats.regNums(regresults)
GetStats.dataNums(regresults, cgnmons)
GetStats.regressSlopes(regresults, MODSUBregresults)
#%%
semidata = ProcessData.getSemivariograms(regresults, MODregresults, MODSUBregresults,   
                                         SurfEst, MODSurfEst, MODSUBSurfEst, 
                                         glats, glons)
#%% prepare the dataframe for making a figure or two below
df = GetData.makeDataframe(glons, glats, regresults, pretot, ncfiles)
GetStats.MaritimeSummary(df)

#%% get geopotential height dataset and create stacks
hgts, hgtlat, hgtlon = GetData.getGeopotentialHeight(ncfiles)

#the areas to make stacks for , N. Am. Northwest and Central Europe
stacksNW, stacksNW_std = ProcessData.makeStacks(np.array([[42.5, 57.5], [-125, -95]]), glons, glats, cgnmons, 
                                                cgtcov, cgseas)
stacksEur, stacksEur_std = ProcessData.makeStacks(np.array([[47.5, 57.5], [5, 20]]), glons, glats, cgnmons, 
                                                  cgtcov, cgseas)

#%% Make figures for manuscript

#figure 1: the trends map 
MakeFigures.changemap(glons, glats, regresults, cgseas, cgtcov, cgnmons, 
                      years, aggplots)

#figure 2: North America and Europe plot
MakeFigures.stackCompositeFig(hgts, hgtlat, hgtlon, years, stacksNW, stacksNW_std, stacksEur, 
                      stacksEur_std, aggplots)


#Figure 3, two panel plot of pressure/omega slope/intercept
MakeFigures.maritime(df, aggplots)


#figure 4, the places where model change matches residence time changes
#this is now the slope of d18O estimate based on residence time changes (slope/slope comparisons)
MakeFigures.PlotModRT(aggplots, Lons, Lats,  
                     -RTMODregresults[:, :, 2, :], MODregresults[:, :, 2, :])

