# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 10:23:07 2019

@author: u0929173
"""
import os
import numpy as np
import netCDF4 as nc4
import copy
import csv
import pandas as pd
import datetime as dt

def getCAM(filename, CAMdir):
    os.chdir(CAMdir)
    
    f = nc4.Dataset(filename, format = 'NETCDF4')
    Lons = f.variables['lon'][:]
    Lats = f.variables['lat'][:]
    
    #these are in m/s
    R18O = np.divide(f.variables['PRECT_H218O'][:], f.variables['PRECT_H216O'][:])
    
    #this is originally in m/s
    Pre = f.variables['PRECT_H2O'][:]*3600*24*1000 #this should now be in kg/m2/day now
    Tmp = f.variables['TREFHT'][:]
    PW = f.variables['TMQ'][:] #kg/m2
    f.close()
    
    d18O = (R18O-1)*1000
    d18O = np.ma.array(d18O, mask = abs(d18O) > 200)
    Pre = np.ma.array(Pre, mask = Pre <= 0.001)
    dates = np.arange(np.datetime64('1940-01'), np.datetime64('2015-01'), 1)
    
    mons = 1+dates.astype(float)%12
    years = dates.astype('datetime64[Y]').astype(float)+1970
    #turn d18O into a z-score
            
    grids = np.zeros((len(Lats), len(Lons), len(np.unique(years)), 2))
    grids.fill(np.nan)
    pregrids = copy.copy(grids)
    tmpgrids = copy.copy(grids)
    pwgrids = copy.copy(grids)
    #index zero is winter, index 1 is summer
    midind = [1, 7]
    for i in range(len(Lats)):
        for j in range(len(Lons)):
            for k in range(2):
                subinds = np.where(mons == midind[k])[0]
                for l in subinds:
                    inds = np.array([l-1, l, l+1])
                    inds = inds[(inds>=0)&(inds<len(mons))]
                    if any(d18O[inds, i, j]):
                        (grids[i, j, int(years[l]-years[0]), k], 
                               pregrids[i, j, int(years[l]-years[0]), k]) = np.average(d18O[inds, i, j], 
                                                                                weights = Pre[inds, i, j], 
                                                                        returned = True)
                        tmpgrids[i, j, int(years[l]-years[0]), k] = np.average(Tmp[inds, i, j])
                        pwgrids[i, j, int(years[l]-years[0]), k] = np.average(PW[inds, i, j])
    flipind = np.where(Lons > 180)[0][0]
    Lons = np.concatenate((Lons[flipind:]-360, Lons[:flipind]))
    d18O = np.concatenate((d18O[:, :, flipind:], d18O[:, :, :flipind]), axis = 2)
    Pre = np.concatenate((Pre[:, :, flipind:], Pre[:, :, :flipind]), axis = 2)
    Tmp = np.concatenate((Tmp[:, :, flipind:], Tmp[:, :, :flipind]), axis = 2)
    PW = np.concatenate((PW[:, :, flipind:], PW[:, :, :flipind]), axis = 2)
    
    grids = np.concatenate((grids[:, flipind:, :, :], grids[:, :flipind, :, :]), axis = 1)
    pregrids = np.concatenate((pregrids[:, flipind:, :, :], pregrids[:, :flipind, :, :]), axis = 1)
    tmpgrids = np.concatenate((tmpgrids[:, flipind:, :, :], tmpgrids[:, :flipind, :, :]), axis = 1)
    pwgrids = np.concatenate((pwgrids[:, flipind:, :, :], pwgrids[:, :flipind, :, :]), axis = 1)
    
    return Lons, Lats, mons, years, dates, d18O, Pre, Tmp, PW, grids, pregrids, tmpgrids, pwgrids

def getdata(ncfiles):

    print('Reading NetCDF4 file...')
    os.chdir(ncfiles)
    f = nc4.Dataset('Precipitation_WI.nc', 'r', format = 'NETCDF4')
    start = copy.copy(f.variables['start'][:])
    start = np.datetime64(start[0])
    t = copy.copy(f.variables['time'][:])
    t = start + t.astype('timedelta64[D]')
    
    iso = copy.copy(f.variables['isoarray'][:])
    Pre = copy.copy(f.variables['precip'][:])
    Tmp = copy.copy(f.variables['temper'][:])
    Lats = copy.copy(f.variables['latitudes'][:])
    Lons = copy.copy(f.variables['longitudes'][:])
    SiteID = copy.copy(f.variables['siteid'][:])
    SiteID = SiteID.astype('|S50')         
    ProjectID = copy.copy(f.variables['projects'][:])
    ProjectID = ProjectID.astype('|S6')
    
    f.close()
    iso[iso <= -9000] = np.nan
    Pre[Pre == -1] = np.nan
    Tmp[Tmp <= -9000] = np.nan
    #these are the bounds of the lat/lon boxes
    boxsize = 5
    glats_edges = np.arange(-90, 95, boxsize)
    glons_edges = np.arange(-177.5, 187.5, boxsize)
    glats = np.average(np.concatenate((glats_edges[1:, None], glats_edges[:-1, None]), axis = 1), axis = 1)
    glons = np.average(np.concatenate((glons_edges[1:, None], glons_edges[:-1, None]), axis = 1), axis = 1)
    #store a latlon index tuple in this array
    latlonind = np.zeros((len(iso), 2))
    
    #get the indices where the isotope records will be placed
    for i in range(len(iso)):
        iy = np.argmin(abs(Lats[i]-glats))
        ix = np.argmin(abs(Lons[i]-glons))
        if (min(abs(Lats[i]-glats)) > boxsize) or (min(abs(Lons[i]-glons)) > boxsize):
            print(Lats[i], Lons[i])
        latlonind[i, :] = np.array([iy, ix])    

    #get annual grids 
    Hfilename = 'Hma.asc'
    Ofilename = 'Oma.asc'        
    j = 0
    #create the grid for matching
    with open(Hfilename, 'r') as f:
        freader = csv.reader(f, delimiter = '\t')
        for row in freader:
            if  j == 0:
                meta = np.array(row)[:, None] 
            elif j <= 5:
                meta = np.concatenate((meta, np.array(row)[:, None]), axis = 1)
            else:
                f.close()
                break
            j = j+1
    

    ncols = int(meta[1, 0])
    nrows = int(meta[1, 1])
    xllcorner = float(meta[1,2])
    yllcorner = float(meta[1, 3])
    cellsize = float(meta[1, 4])
    fillval = float(meta[1, 5])
            
    gridlats_e = (np.arange(nrows, -1, -1))*cellsize+yllcorner
    gridlons_e = (np.arange(ncols))*cellsize+xllcorner
    gridlats = np.average(np.concatenate((gridlats_e[1:, None], gridlats_e[:-1, None]), axis = 1), axis = 1)
    gridlons = np.average(np.concatenate((gridlons_e[1:, None], gridlons_e[:-1, None]), axis = 1), axis = 1)
    #make the array to fill

    isogrids = np.zeros((1, nrows, ncols, 2))
    with open(Hfilename, 'r') as f:  
        freader = csv.reader(f, delimiter = '\t')
        for _ in range(7):
            row = next(freader)                             
        isogrids[0, :, :, 0]  = np.reshape(np.array(row[:-1]).astype(float), (nrows, ncols), order = 'C')
    f.close()
   #oxygen is index 1
    with open(Ofilename, 'r') as f:  
        j = 0
        freader = csv.reader(f, delimiter = '\t')
        for _ in range(7):
            row = next(freader)            
        isogrids[0, :, :, 1]  = np.reshape(np.array(row[:-1]).astype(float), (nrows, ncols), order = 'C')
    f.close()
    isogrids[isogrids == fillval] = np.nan        
    
    
    coarsegridavg = np.zeros((len(glons), len(glats), 2))
    coarsegridavg.fill(np.nan)
    for i in range(len(glats)):
        #no need to fill anything above the maximum latitude
        if glats[i] <= np.max(gridlats):
            iy = np.where(((glats[i]+boxsize/2.0)>= gridlats)&((glats[i]-boxsize/2.0) < gridlats))[0]
            for j in range(len(glons)):
                ix = np.where(((glons[j]+boxsize/2.0) >= gridlons)&((glons[j]-boxsize/2.0) < gridlons))[0]
                vals = isogrids[0, min(iy):max(iy), min(ix):max(ix), :]
                if vals.all() != -9999:                    
                    vals = np.where(vals == -9999, np.nan, vals)
                    coarsegridavg[j, i, :] = np.nanmean(np.nanmean(vals, axis = 0), axis = 0)
            
    #calculate anomalies based on the closest grid value
    #hydrogen is index 1 and oxygen is index 2
    isoanom = np.zeros(iso.shape)
    #go by site
    seasvals = np.zeros((iso.shape[0], 2))
    seasvals.fill(np.nan)
    monthlist = t.astype('datetime64[M]').astype(float)%12
    searchsize = 15
    canaggflag = np.zeros(iso.shape[0])
    for i in range(iso.shape[0]):
        #find the isogrid cell that matches the site location    
        iy = np.argmin(abs(gridlats - Lats[i]))
        ix = np.argmin(abs(gridlons - Lons[i]))  
        seasvals[i, :] = isogrids[:, iy, ix, :] 
        if np.isnan(seasvals[i]).any():
            
            #make a mask of the available areas
            mask = np.isnan(isogrids[0, (iy-searchsize):(iy+searchsize), (ix-searchsize):(ix+searchsize), 1])
            #calculate euclidian distance to each grid cell
            y, x = np.meshgrid(gridlats[(iy-searchsize):(iy+searchsize)], gridlons[(ix-searchsize):(ix+searchsize)])
            dist = np.sqrt(np.add(np.power((x-Lons[i]), 2), np.power((y-Lats[i]), 2)))
            dist = np.ma.array(dist, mask = mask)
            xind, yind = np.where(dist == np.min(dist))
            minx = x[:, 0][xind[0]]
            miny = y[:, 0][yind[0]]
            
            iy = np.argmin(abs(gridlats - miny))
            ix = np.argmin(abs(gridlons - minx))          
            seasvals[i, :] = isogrids[:, iy, ix, :]            
    
        for j in range(12):
            inds = np.where(j == monthlist)[0]
            if np.isnan(seasvals[i]).any():
                #index 1 means that it cannot be aggregated to other sites
                isoanom[i, inds, 1:3] = iso[i, inds, 1:3]
                canaggflag[i] = 1
            else:
                #index zero means it can be aggregated to other sites
                isoanom[i, inds, 1:3] = np.subtract(iso[i, inds, 1:3], seasvals[i, :])
                #print(isoanom[i, inds, 1:3][~np.isnan(isoanom[i, inds, 1:3])])
                canaggflag[i] = 0
    
    coarsegrid = np.zeros((len(glats), len(glons), len(t), 2))
    projgrid = np.empty((len(glats), len(glons), len(t)), dtype='S50')
    pregrid = np.zeros((len(glats), len(glons), len(t)))
    tmpgrid = copy.copy(pregrid)
    gridDelT = np.zeros((len(glats), len(glons)))
    ninds = np.zeros((len(glats), len(glons)))
    coarsegrid.fill(np.nan)
    projgrid.fill('')
    indlist = np.array([], dtype='int')
    for i in range(coarsegrid.shape[0]):
        for j in range(coarsegrid.shape[1]):
            inds = np.where((latlonind[:, 1] == j) & (latlonind[:, 0] == i))[0]
            #need to remove the non-aggregatable sites
            if any(inds):           
                if (len(inds) > 1) and sum(canaggflag[inds]) > 0 and sum(canaggflag[inds])<len(inds):
                    print(len(inds))
                    print('removing non-aggregatable sites {}'.format(inds[canaggflag[inds]==1]))
                    inds = inds[canaggflag[inds]==0]
                    print(len(inds))
                #print(len(inds))
                #now aggregate the sites
                #may also want to make a diagnostic plot of the sites that are part of the integration
                if len(inds) > 1:
                    #add back the grid scale center value
                    coarsegrid[i, j, :, :] = np.nanmean(np.add(isoanom[inds, :, 1:3], coarsegridavg[j, i, :]), axis = 0)
                    pregrid[i, j, :] = np.nanmean(Pre[inds, :], axis = 0)
                    tmpgrid[i, j, :] = np.nanmean(Tmp[inds, :], axis = 0)
                    for k in range(len(t)):
                        if any(ProjectID[inds, k]):
                            prj = list(np.unique(ProjectID[inds, k][ProjectID[inds, k] != b'']))
                            print(prj)
                            if len(prj) == 1:
                                print(prj[0])
                                projgrid[i, j, k] = prj[0]
                            else:
                                prjs = b', '.join(prj)
                                print(prjs)
                                projgrid[i, j, k] = prjs
                #this means that len(inds) == 1
                else:
                    if (~np.isnan(coarsegridavg[j, i, :]).all()) and (canaggflag[inds] == 0):
                        coarsegrid[i, j, :, :] = np.add(isoanom[inds, :, 1:3], coarsegridavg[j, i, :])  
                    else:
                        coarsegrid[i, j, :, :] = isoanom[inds, :, 1:3]
                    pregrid[i, j, :] = Pre[inds, :]
                    tmpgrid[i, j, :] = Tmp[inds, :]
                    for k in range(len(t)):
                        projgrid[i, j, k] = ProjectID[inds[0], k]
                
                if ~np.isnan(isoanom[inds,:,1:3]).any(): 
                    mask = ~np.isnan(coarsegrid[i, j, :, 1])
                    gridDelT[i, j] = np.max(t.astype('datetime64[M]').astype(float)[mask]) - np.min(t.astype('datetime64[M]').astype(float)[mask])
                #return a list of the station indices
                indlist = np.concatenate((indlist, inds))
    return t, glons, glats, coarsegrid, pregrid, gridDelT, ninds, tmpgrid, indlist, projgrid   

def getseasonalgrids(glats, glons, coarsegrid, pregrid, t, tmpgrid, indlist, projgrid):
    #this part gets the seasonal average value for each coarse grid cell, then gets the regression slope
    seasdict = {'DJF': 1, 'JJA': 7}
    cgseas = np.zeros((len(glats), len(glons), len(np.unique(t.astype('datetime64[Y]').astype(float))), 2))
    cgseas.fill(np.nan)
    cgpre = np.zeros((len(glats), len(glons), len(np.unique(t.astype('datetime64[Y]').astype(float))), 2))
    pretot = np.zeros((len(glats), len(glons), len(np.unique(t.astype('datetime64[Y]').astype(float))), 2))
    cgtmp = copy.copy(cgpre)
    cgnmons = np.zeros((len(glats), len(glons), 2))
    cgtcov = np.zeros((len(glats), len(glons), 2))
    cgprjs = np.empty((len(glats), len(glons), len(np.unique(t.astype('datetime64[Y]').astype(float))), 2), dtype='S50')
    cgprjs.fill('')    
    years = t.astype('datetime64[Y]').astype(float)+1970
    mons = t.astype('datetime64[M]').astype(float)%12
    for season, seasind in zip(seasdict.keys(), range(len(seasdict.keys()))):
        monmatch = seasdict[season]
        for (y, i) in zip(np.unique(years), range(len(np.unique(years)))):
            subind = np.where((years == y)&(mons == monmatch))[0]
            subinds = np.array([subind-1, subind, subind+1])
            subinds = subinds[(subinds>=0)&(subinds < len(years))]
            #this is getting just oxygen
            vals = np.ma.array(coarsegrid[:, :, subinds, 1], mask = np.isnan(coarsegrid[:, :, subinds, 1]))
            
            wts = np.ma.array(pregrid[:, :, subinds], mask = np.isnan(coarsegrid[:, :, subinds, 1]))
            pretot[:, :, i, seasind] = np.sum(pregrid[:, :, subinds], axis = 2)
            cgtmp[:, :, i, seasind] = np.average(tmpgrid[:, :, subinds], axis = 2)
            cgseas[:, :, i, seasind], cgpre[:, :, i, seasind] = np.ma.average(vals, weights = wts, axis = 2, returned = True) 
            #apply the precipitation threshold
            cgseas[:, :, i, seasind] = np.where(np.divide(cgpre[:, :, i, seasind], pretot[:, :, i, seasind]) >= 0.5, cgseas[:, :, i, seasind], np.nan)
            for loni in range(len(glons)):
                for lati in range(len(glats)):
                    if any(projgrid[lati, loni, subinds]):
                        prj = list(np.unique(projgrid[lati, loni, subinds][projgrid[lati, loni, subinds] != b'']))
                        if len(prj) ==1:
                            cgprjs[lati, loni, i, seasind] =  prj[0]
                        else:
                            prjs = b', '.join(prj)
                            prjs = np.unique(prjs)
                            cgprjs[lati, loni, i, seasind] = prjs[0]
            
    #apply thresholds to cgseas
    cgnmons = np.nansum(~np.isnan(cgseas), axis = 2)
    for i in range(cgtcov.shape[0]):
        for j in range(cgtcov.shape[1]):
            for k in range(2):
                mask = ~np.isnan(cgseas[i, j, :, k])
                if any(mask):
                    cgtcov[i, j, k] = np.max(np.unique(years)[mask])-np.min(np.unique(years)[mask])
    #note that cgnseas is cgnmons in the rest of the code. 
    return cgseas, cgpre, pretot, cgnmons, cgtcov, years, mons, cgtmp, cgprjs

def makeDataframe(glons, glats, regresults, pretot, ncfiles):
    #locations of things
  
    pres = nc4.Dataset(os.path.join(ncfiles, "regrid_Pres_surf_1960_2014.nc"))
    nclons = pres.variables['lon'][:]
    nclats = pres.variables['lat'][:]
    pres.close()
    
    columns = ['Lat', 'Lon', 'Season', 'Slope', 'Slope_std', 'Maritime', 'Clim', 'Press', 
           'Omega', 'Omega_std', 'Press_beta', 'Press_std', 'T', 'T_std', 'P']
    
    inds = np.where(~np.isnan(regresults[:, :, 2, :]))
    df = pd.DataFrame(columns = columns)
    precip = np.nanmean(pretot, axis = 2)
    os.chdir(ncfiles)
    omegareg = np.load("omegareg_re.npy")
    presreg = np.load("presreg_re.npy")
    treg = np.load("treg.npy")
    
    for row, i, j, k in zip(range(len(inds[0])), inds[0], inds[1], inds[2]):
        df.loc[row, 'Lat'] = glats[i]
        df.loc[row, 'Lon'] = glons[j]
        if k == 0:
            df.loc[row, 'Season'] = 'DJF'
        else:
            df.loc[row, 'Season'] = 'JJA'
            
        df.loc[row, 'Slope'] = regresults[i, j, 2, k]
        df.loc[row, 'Slope_std'] = regresults[i, j, 3, k]
        df.loc[row, 'P'] = precip[i, j, k]
        
        if glons[j] < 0:
            lonval = glons[j]+360
        loni = np.where(lonval == nclons)[0][0]
        lati = np.where(glats[i] == nclats)[0][0]    
        
        df.loc[row, 'Press'] = presreg[i, j, 2, k]
        df.loc[row, 'Press_std'] = presreg[i, j, 3, k]
        df.loc[row, 'Omega'] = omegareg[i, j, 2, k]
        df.loc[row, 'Omega_std'] = omegareg[i, j, 3, k] 
        df.loc[row, 'T'] = treg[lati, loni, 2, k]
        df.loc[row, 'T_std'] = treg[lati, loni, 3, k]
        
    #koppen classification key:
    #http://koeppen-geiger.vu-wien.ac.at/shifts.htm
    #Main Climates:
    #   A: Equatorial, B: Arid C: warm temperate D: snow E: Polar
    #Precipitation
    # W: desert S: Steppe f:fully humid s:summer dry w:winter dry m:monsoonal
    #Temperature
    #   h:hot arid k:cold arid a:hot summer b:warm summer c:cool summer 
    #   d:extremely continental F: polar frost T: polar tundra
    
    koppen = pd.read_csv('1976-2000_ASCII.txt', header = 0, delim_whitespace=True)
    res = 2.5
    maxgrids = 100
    koppenregrid = pd.DataFrame(columns = ['Lat', 'Lon', 'ngrids', 'Cls'])
    i = 0
    for y in glats:
        for x in glons:
            koppenregrid.loc[i, 'Lat'] = y
            koppenregrid.loc[i, 'Lon'] = x
            
            mask1 = (koppen['Lat'] > (y - res))&(koppen['Lat'] < (y + res))
            mask2 = (koppen['Lon'] > (x - res))&(koppen['Lon'] < (x + res))
            mask = mask1&mask2
            ngrids = sum(mask)
            koppenregrid.loc[i, 'ngrids'] = ngrids
            
            if ngrids > 0:
                classes, count  = np.unique(koppen['Cls'].values[mask], return_counts = True)
                koppenregrid.loc[i, 'Cls'] = classes[np.argmax(count)]            
                
            i = i+1
            
    for i in range(len(df)):    
        koppenind = np.where((df.loc[i, 'Lon'] == koppenregrid['Lon'].values)&(df.loc[i, 'Lat'] == koppenregrid['Lat'].values))[0][0]    
        #if some or most of grid is in ocean, we will consider it maritime
        #need to mask out the mediterranian area
        if ((((koppenregrid.loc[koppenind, 'Lat'] > 30)&(koppenregrid.loc[koppenind, 'ngrids'] < 70))|
        ((koppenregrid.loc[koppenind, 'Lat'] <= 30)&(koppenregrid.loc[koppenind, 'ngrids'] < 90)))&
         ~((koppenregrid.loc[koppenind, 'Lat'] > 30)&((koppenregrid.loc[koppenind, 'Lon'] > 9)&(koppenregrid.loc[koppenind, 'Lon'] < 57)))):     
            df.loc[i, 'Maritime'] = 1
        else:
            df.loc[i, 'Maritime'] = 0
        if koppenregrid.loc[koppenind, 'ngrids'] > 0:
            df.loc[i, 'Clim'] = koppenregrid.loc[koppenind, 'Cls']
        else:
            df.loc[i, 'Clim'] = 'N'
            
    return df

def getGeopotentialHeight(ncfiles):
    #this only gets winter values
    f = nc4.Dataset(os.path.join(ncfiles, 'hgt500mon.nc'), 'r', format = 'NETCDF4')
    timedel = f['time'][:]
    basetime = dt.datetime(1800, 1, 1)
    hgttime = [basetime+dt.timedelta(hours = th) for th in timedel]
    hgtlon = f['lon'][:]
    hgtlat = f['lat'][:]
    hgt = f['hgt'][:]
    f.close()
    
    #collapse to only winter values between 1960 and 2016
    mons = [th.month for th in hgttime]
    inds = np.where(np.array(mons) == 12)[0]
    inds = inds[inds < len(mons)-2]
    hgts = np.zeros((len(inds), len(hgtlat), len(hgtlon)))
    for i, ind in zip(range(len(inds)), inds):
        hgts[i, :, :] = np.average(hgt[[ind, ind+1, ind+2], :, :], axis = 0)
    return hgts, hgtlat, hgtlon

def getTeleIndex(teleloc, telename):
    INDEX = pd.read_csv(os.path.join(teleloc, telename), sep = ' ', header = 0)
    INDEX = INDEX[(INDEX.loc[:, 'Year'] >=1959)&(INDEX.loc[:, 'Year'] <=2016)].reset_index()
    end = len(INDEX)-2
    INDEX_DJF = np.average(np.concatenate((INDEX.loc[:end,['Dec']].values, 
                                         INDEX.loc[1:, ['Jan', 'Feb']].values), 
                            axis = 1), axis = 1)
    return INDEX_DJF
