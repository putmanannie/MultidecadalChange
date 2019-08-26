# -*- coding: utf-8 -*-
"""
Script that does all processing, analysis, makes figures, and produces statistics 
for the paper. 
"""
import os
import statsmodels.api as sm
import scipy as sp
import numpy as np
import netCDF4 as nc4
import copy
import csv
import matplotlib as mpl
import pandas as pd
import matplotlib.patches as patches
import cartopy.crs as ccrs
import cartopy as ctp
import datetime as dt
import urllib

rootdir = """C:\Users\u0929173\Documents\Python\WaterIsotopeAnalysis_012016\GeospatialInterpolation\Manuscript"""
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
CAMdir = r"""W:\iCAM\forAnnie"""  #this is big so needs to be held elsewhere
getregrid = False
#%% ######## GET AND PROCESS DATA (MODEL AND MEASUREMENT) ##################

#Get and process model data (3 ensemble members)
(Lons, Lats, icyears, d18O, Pre, Tmp, PW, grids, 
     pregrids, tmpgrids, pwgrids) = ProcessData.processCAM(os, np, nc4, copy, GetData, CAMdir)

#Get and process measurement data
(glats, glons, cgseas, cgpre, pretot, cgnmons, 
     cgtcov, years, mons) = ProcessData.processData(os, np, nc4, csv, copy, GetData)

startind = np.where(np.unique(icyears) == np.min(years))[0][0]
endind = np.where(np.unique(years) == np.max(icyears))[0][0]+1
gridsub = np.where(np.isnan(cgseas[:, :, :endind, :]), cgseas[:, :, :endind, :], grids[:, :, startind:, :])
MODSUBSurfEst, MODSUBSurfEstSE = EstimateRT.getSurfEst(np, gridsub)
#use the same time frame
MODSurfEst, MODSurfEstSE = EstimateRT.getSurfEst(np, grids[:, :, startind:, :])
#%% ######### RESIDENCE TIME ESTIMATE OF D18O TIMESERIES (MODEL AND MEASUREMENT) ################

#estimate residence time d18O timeseries from model data
d18O_ts, yearrange = EstimateRT.EstModRT(np, copy, icyears, grids, pregrids, pwgrids)

#estimate residence time for measurement from daily scale reanalysis data
SurfEst, SurfEstSE = EstimateRT.getSurfEst(np, cgseas)
if getregrid:
    EstimateRT.GetRTData(os, np, dt, urllib, prefileindir)
    EstimateRT.regridRTdata(iris, np, os)
#this is returning a lot of nans... something is wrong
Tseas, d18O_est = EstimateRT.getMeasRT(os, np, nc4, copy, SurfEst, SurfEstSE, prefileindir, ncfiles)
cgseasRT = np.swapaxes(np.swapaxes(d18O_est, 2, 0), 3, 1)

#%% ############# DO REGRESSIONS OF MODEL, DATA, AND RESIDENCE TIME TIMESERIES ########
#These should come out in per mil per decade
#data
regresults = GetChange.GetMultidecadalChange(np, sm, glats, 
                                             glons, years, cgseas, 
                                             cgtcov, cgnmons)
#model d18O
MODregresults = GetChange.getCAMChange(np, sm, icyears, years, Lats, 
                             Lons, grids, cgnmons, cgtcov, [0, 0])
#subset
MODSUBregresults = GetChange.GetMultidecadalChange(np, sm, glats, 
                                             glons, years[years <= np.max(icyears)],  
                                             gridsub, cgtcov, cgnmons)

#model residence time
RTMODregresults = GetChange.getCAMChange(np, sm, yearrange, years, Lats, Lons, 
                                         d18O_ts, cgtcov, cgnmons, [0,0])

#reanalysis residence time
RTregresults = GetChange.GetMultidecadalChange(np, sm, glats, 
                                             glons, years, cgseasRT, 
                                             cgtcov, cgnmons)
#summarize comparison of residence time estimate and model or measurement slope estimate
GetStats.summary(np, RTregresults, regresults)
GetStats.summary(np, RTMODregresults, MODregresults)
#%%
semidata = ProcessData.getSemivariograms(os, np, mpl, sp, pd, copy, ProcessData.distance, 
                                         ProcessData.fitSemiVariogram, regresults, 
                                         MODregresults, MODSUBregresults, SurfEst,  
                                         MODSurfEst, MODSUBSurfEst, glats, glons)
#%% prepare the dataframe for making a figure or two below
df = GetData.makeDataframe(os, np, pd, nc4, glons, glats, regresults, pretot, ncfiles)
GetStats.MaritimeSummary(np, sm, df)

#%% get geopotential height dataset and create stacks
hgts, hgtlat, hgtlon = GetData.getGeopotentialHeight(os, np, dt, nc4, ncfiles)
PNA_DJF = GetData.getTeleIndex(os, pd, np, teleloc, 'PNA.txt')
NAO_DJF = GetData.getTeleIndex(os, pd, np, teleloc, 'NAO.txt')
#the areas to make stacks for , N. Am. Northwest, Southeast and Central Europe
NWarea = np.array([[42.5, 57.5], [-125, -95]])
SEarea = np.array([[32.5, 52.5], [-95, -60]])
Eur = np.array([[47.5, 57.5], [5, 20]])
stacksNW, stacksNW_std = ProcessData.makeStacks(np, NWarea, glons, glats, cgnmons, 
                                                cgtcov, cgseas)
stacksSE, stacksSE_std = ProcessData.makeStacks(np, SEarea, glons, glats, cgnmons, 
                                                cgtcov, cgseas)
stacksEur, stacksEur_std = ProcessData.makeStacks(np, Eur, glons, glats, cgnmons, 
                                                  cgtcov, cgseas)

#%%
#figure 1: the change map
MakeFigures.changemap(os, mpl, np, patches, ctp, ccrs, glons, glats, 
                   regresults, cgseas, cgtcov, cgnmons, years, aggplots)
#figure 2 or Figure S1, the places where model change matches residence time changes
#this is now the slope of d18O estimate based on residence time changes (slope/slope comparisons)
MakeFigures.PlotModRT(os, np, mpl, ctp, ccrs, aggplots, Lons, Lats,  
                     RTMODregresults[:, :, 2, :], MODregresults[:, :, 2, :])
#Figure 2, stacked bar plot
MakeFigures.StackedBar(os, np, mpl, patches, df, aggplots)

#Figure of semivariogram lag distances
MakeFigures.MakeSemivariogramFig(os, np, mpl, semidata, aggplots)

#Figure 3, two panel plot of pressure/omega slope/intercept
MakeFigures.maritime(os, sm, mpl, df, aggplots)

#make geopotential height figures
MakeFigures.MakeGeopotentialRegress(os, np, mpl, ctp, ccrs, stacksNW, hgts, 
                                    hgtlat, hgtlon, 'NW', '500')
MakeFigures.MakeGeopotentialRegress(os, np, mpl, ctp, ccrs, stacksEur, hgts, 
                                    hgtlat, hgtlon, 'Eur', '500')
MakeFigures.MakeStackFig(os, np, mpl, years, stacksNW, stacksNW_std, stacksSE, 
                         stacksSE_std, PNA_DJF, aggplots)
