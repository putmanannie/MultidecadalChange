# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 12:08:25 2019

@author: u0929173
"""
import numpy as np
import copy
import datetime as dt
import urllib
import os
#import iris
import netCDF4 as nc4

def EstModRT(icyears, grids, pregrids, pwgrids):
    #these are the long term values bookending the period of reference.
    startyr = 1960
    endyr = 2015
    yearrange = np.arange(startyr, endyr)
    d18O_ts = np.zeros((grids.shape[0], grids.shape[1], len(yearrange), 2))
    Dd18O_ts = np.zeros((grids.shape[0], grids.shape[1], len(yearrange), 2))
    
    #get the long term avg values
    sind = np.argmin(abs(np.unique(icyears)-startyr))
    eind = np.argmin(abs(np.unique(icyears)-endyr))
    premean = np.nanmean(pregrids[:, :, sind:eind, :], axis = 2) 
    pwmean = np.nanmean(pwgrids[:, :, sind:eind, :], axis = 2)
    gridmean= np.nanmean(grids[:, :, sind:eind, :], axis = 2)
    
    for (i, y) in zip(range(len(yearrange)), yearrange):
        subind = np.where(np.unique(icyears) == y)[0][0]
        Qanom = np.divide(np.subtract(pwgrids[:, :, subind, :], pwmean), pwmean)
        Panom = np.divide(np.subtract(pregrids[:, :, subind, :], premean), premean)
        d18O_ts[:, :, i, :] = np.add(np.multiply(-0.4*np.subtract(Qanom, Panom), gridmean), gridmean)
        Dd18O_ts[:, :, i, :] = -0.4*np.subtract(Qanom, Panom)
    return d18O_ts, Dd18O_ts, yearrange

def GetRTData(prefileindir):
    #get the data if it needs to be gotten, probs wont after this time
    getncfiles = False
    
    if getncfiles:
        thisyear = dt.datetime.now().year
        prwtfilenames = ["pr_wtr.eatm."+str(year)+".nc" for year in range(1960, thisyear+1)]
        pratefilenames = ["prate.sfc.gauss."+str(year)+".nc" for year in range(1960, thisyear+1)]
        tmpfilenames = ["air.sig995."+str(year)+".nc" for year in range(1960, thisyear+1)]
        #check if all files are downloaded, if not, ftp to most recent

        for (pwfn, prfn, tpfn) in zip(prwtfilenames, pratefilenames, tmpfilenames):
            if pwfn not in prefileindir:
                urllib.urlretrieve("""ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface/{}""".format(pwfn),
                                   os.path.join(prefileindir, pwfn))
                urllib.urlcleanup()
                
            if prfn not in prefileindir:
                urllib.urlretrieve("""ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface_gauss/{}""".format(prfn),
                                   os.path.join(prefileindir, prfn))
                urllib.urlcleanup()
                
            if tpfn not in prefileindir:
                urllib.urlretrieve("""ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface/{}""".format(tpfn),
                                   os.path.join(prefileindir, tpfn))
                urllib.urlcleanup()
 
#needs iris to be imported but iris is giving issues with this version of python               
def regridRTdata():
    prefiledir = """W:\Annie\Reanalysis"""
    regridfiledir = """W:\Annie\Reanalysis\Regridded"""
    prefileindir = np.sort(os.listdir(prefiledir))
    
    #make the grid to regrid to
    filename = """f.e13.FiPIPDC5CN.5x5.ens1.cam.h0.194001-201412.nc"""
    gridfiledir = r"""W:\iCAM\forAnnie"""
    getgrid1 = iris.load(os.path.join(gridfiledir, filename))
    basegrid = getgrid1[0]
    basegrid.coord('longitude').guess_bounds()
    basegrid.coord('latitude').guess_bounds()
    #mdtol = percent of masked data in target grid box. Currently 50% can be missing
    #anything less than 1 means that this will not be mass conservative, fyi
    scheme = iris.analysis.AreaWeighted(mdtol = 0.5)
    
    for f in prefileindir:
        if f[:5] == 'air.s':
            varname = u'mean Daily Air temperature at sigma level 995'
            cube = iris.load_cube(os.path.join(prefiledir, f), varname)
        elif f[:5] == 'pr_wt':
            varname = u'mean Daily Precipitable Water for entire atmosphere'
            cube = iris.load_cube(os.path.join(prefiledir, f), varname)        
        else:
            varname = u'mean Daily Precipitation Rate at surface'
            cube = iris.load_cube(os.path.join(prefiledir, f), varname)
            
        cube.coord('longitude').guess_bounds()
        cube.coord('latitude').guess_bounds()
        regridded_cube = cube.regrid(basegrid, scheme) 
        iris.save(regridded_cube, os.path.join(regridfiledir, 'regrid_{}'.format(f)))

def getSurfEst(cgseas):
    SurfEst = np.nanmean(cgseas, axis = 2)
    sqrt_nyrs = np.sqrt(np.sum(~np.isnan(cgseas), axis = 2))
    sqrt_nyrs = np.where(sqrt_nyrs == 0, np.nan, sqrt_nyrs)
    SurfEstSE = np.divide(np.nanstd(cgseas, axis = 2), sqrt_nyrs)
    return SurfEst, SurfEstSE

def getMeasRT(SurfEst, SurfEstSE, prefileindir, ncfiles): 
    #indices corresponding to july and january.
    monind = [0, 6]
    startyr = 1960
    endyr = 2016
    yearrange = np.arange(startyr, endyr+1)

    initialize = True

    time = np.array([], dtype = 'datetime64[D]')
    #make the whole record for the timeperiod then mask it/average it at the end
    for y in range(startyr, endyr+1):       
        
        #open the precipitation rate file
        dspr = nc4.Dataset(os.path.join(prefileindir, 'Regridded', 'regrid_prate.sfc.gauss.{}.nc'.format(y)))
        #get precipitable water
        dspw = nc4.Dataset(os.path.join(prefileindir, 'Regridded', 'regrid_pr_wtr.eatm.{}.nc'.format(y)))
        tmpr = nc4.Dataset(os.path.join(prefileindir, 'Regridded', 'regrid_air.sig995.{}.nc'.format(y)))   
        
        if initialize:
            
            lat = dspr.variables['lat'][:]
            lon = dspr.variables['lon'][:]

            P = np.ma.empty((len(yearrange), 12, len(lat), len(lon)))
            Q = copy.copy(P)
            T = copy.copy(P)
            initialize = False
            
        #this is how we handle time
        #time is same in all files
        t = dspr.variables['time'][:]/24
        t = t.astype('timedelta64[D]')
        basetime = np.datetime64('1800-01-01', 'D')
        time= t+basetime
        
        #get precipitation rate and close
        #this is in kg/m**2/s (or mm/s)
        prate = dspr.variables['prate'][:]
#        prate = np.where(prate < 0, 0, prate)
        #this is in kg/m**2
        prwt = dspw.variables['pr_wtr'][:]
#        prwt = np.where((prwt < 1000)&(prwt > 0), prwt, np.nan)

        tmp = tmpr.variables['air'][:]
#        tmp = np.where(tmp < 1000, tmp, np.nan)  
        #get monthly values
        for m in range(12):
            ind = y-startyr
            subinds = np.where(m == time.astype('datetime64[M]').astype(int)%12)[0]
            #these are all masked arrays, so use np.ma
            #When assigned to P, a masked array, the mask stays
            P[ind, m, :, :] = np.ma.mean(prate[subinds], axis = 0)
            T[ind, m, :, :] = np.ma.mean(tmp[subinds], axis = 0)
            Q[ind, m, :, :] = np.ma.mean(prwt[subinds], axis = 0)
           
        tmpr.close()            
        dspw.close()
        dspr.close()  
        

    #Need seasonal grids now
    Pseas = np.ma.empty((P.shape[0], 2, P.shape[2], P.shape[3]))
    Qseas = copy.copy(Pseas)
    Tseas = copy.copy(Pseas)
    for y in range(P.shape[0]):
        for s, i in zip(monind, range(len(monind))):            
            Pseas[y, i, :, :] = np.ma.mean(P[y, [(s-1), s, (s+1)], :, :], axis = 0)
            Qseas[y, i, :, :] = np.ma.mean(Q[y, [(s-1), s, (s+1)], :, :], axis = 0)
            Tseas[y, i, :, :] = np.ma.mean(T[y, [(s-1), s, (s+1)], :, :], axis = 0)
###################################################
    #need to flip the values we use
    flipind = np.where(lon > 180)[0][0]
    lon = np.concatenate((lon[flipind:]-360, lon[:flipind]))

    Qseas = np.concatenate((Qseas[:, :, :, flipind:], Qseas[:, :, :, :flipind]), axis = 3)
    Pseas = np.concatenate((Pseas[:, :, :, flipind:], Pseas[:, :, :, :flipind]), axis = 3)
    Tseas = np.concatenate((Tseas[:, :, :, flipind:], Tseas[:, :, :, :flipind]), axis = 3)
    
    Qmean = np.ma.mean(Qseas, axis = 0)
    Pmean = np.ma.mean(Pseas, axis = 0)
    Tmean = np.ma.mean(Tseas, axis = 0)    
      
    d18O_est = np.zeros((len(yearrange), 2, len(lat), len(lon)))
    Dd18O_est = np.zeros((len(yearrange), 2, len(lat), len(lon)))
    for i in (yearrange-yearrange[0]):
            #using subtract and divide with a masked array keeps the mask
            DQ = np.divide(np.subtract(Qseas[i], Qmean), Qmean)
            DP = np.divide(np.subtract(Pseas[i], Pmean), Pmean)
            DRT = np.subtract(DQ, DP)
            #change units of precipitation so it comes out in days
            RT  = np.divide(Qseas[i], Pseas[i]*3600*24)
    
            Dd18O = -0.4*DRT
            Dd18O_est[i, :, :] = Dd18O
            swapSurf = np.swapaxes(np.swapaxes(SurfEst, 2, 0), 2, 1)
            #propagate the standard error of the mean also
            d18O_est[i, :, :] = np.add(np.multiply(Dd18O, swapSurf), swapSurf)
            
    #save the calculated values to a netcdf file
    os.chdir(ncfiles)
    f = nc4.Dataset('RTmodel_d18O.nc', 'w', format = 'NETCDF4')
    ntimes = f.createDimension('times', np.shape(d18O_est)[0])
    lat = f.createDimension('lat', np.shape(Tseas)[2])
    lon = f.createDimension('lon', np.shape(Tseas)[3])
    seasons = f.createDimension('seasons', 2)
    
    TMP = f.createVariable('tmp', 'f', ('times', 'seasons', 'lat', 'lon'))
    TMP[:] = Tseas
    
    d18OEST = f.createVariable('d18Oest', 'f', ('times', 'seasons', 'lat', 'lon'))
    d18OEST[:] = d18O_est
    
    Dd18OEST = f.createVariable('Dd18Oest', 'f', ('times', 'seasons', 'lat', 'lon'))
    Dd18OEST[:] = Dd18O_est
    
    f.close()
    
    return Tseas, d18O_est, Dd18O_est

def loadMeasRT(dataloc):
    f = nc4.Dataset(os.path.join('RTmodel_d18O.nc'), 'r', format = 'NETCDF4')
    Tseas = f.variables['tmp'][:]
    d18O_est = f.variables['d18Oest'][:]
    return Tseas, d18O_est


