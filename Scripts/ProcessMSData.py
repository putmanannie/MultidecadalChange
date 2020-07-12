# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:39:30 2019

@author: u0929173
"""
import numpy as np
import copy
import pandas as pd
import scipy as sp

def processCAM(GetData, CAMdir):
    filename = """f.e13.FiPIPDC5CN.5x5.ctrl.cam.h0.194001-201412.nc"""
    (Lons, Lats, icmons, 
        icyears, dates, d18O0, 
        Pre0, Tmp0, PW0, grids0, 
        pregrids0, tmpgrids0, pwgrids0) = GetData.getCAM(filename, CAMdir)
        
    filename = """f.e13.FiPIPDC5CN.5x5.ens1.cam.h0.194001-201412.nc"""
    (Lons, Lats, icmons, 
        icyears, dates, d18O1, 
        Pre1, Tmp1, PW1, grids1, 
        pregrids1, tmpgrids1, pwgrids1) = GetData.getCAM(filename, CAMdir)
    
    filename = """f.e13.FiPIPDC5CN.5x5.ens2.cam.h0.194001-201412.nc"""
    (Lons, Lats, icmons, 
        icyears, dates, d18O2, 
        Pre2, Tmp2, PW2, grids2, 
        pregrids2, tmpgrids2, pwgrids2) = GetData.getCAM(filename, CAMdir)
    
    d18O =  np.nanmean(np.concatenate((d18O0[None], d18O1[None], d18O2[None]), axis = 0), axis = 0)   
    Pre =  np.nanmean(np.concatenate((Pre0[None], Pre1[None], Pre2[None]), axis = 0), axis = 0) 
    Tmp =  np.nanmean(np.concatenate((Tmp0[None], Tmp1[None], Tmp2[None]), axis = 0), axis = 0) 
    PW =  np.nanmean(np.concatenate((PW0[None], PW1[None], PW2[None]), axis = 0), axis = 0) 
    
    grids = np.nanmean(np.concatenate((grids0[None], grids1[None], grids2[None]), axis = 0), axis = 0)
    pregrids = np.nanmean(np.concatenate((pregrids0[None], pregrids1[None], pregrids2[None]), axis = 0), axis = 0)  
    tmpgrids = np.nanmean(np.concatenate((tmpgrids0[None], tmpgrids1[None], tmpgrids2[None]), axis = 0), axis = 0)  
    pwgrids = np.nanmean(np.concatenate((pwgrids0[None], pwgrids1[None], pwgrids2[None]), axis = 0), axis = 0) 
    
    return Lons, Lats, icyears, d18O, Pre, Tmp, PW, grids, pregrids, tmpgrids, pwgrids

def processData(GetData, ncfiles):
    (t, glons, glats, coarsegrid, 
     pregrid, gridDelT, ninds, tmpgrid, indlist, projgrid) = GetData.getdata(ncfiles)

    (cgseas, cgpre, pretot, 
     cgnmons, cgtcov, years, mons, cgtmp, cgprjs) = GetData.getseasonalgrids(glats, glons, 
                                        coarsegrid, pregrid, t, tmpgrid, indlist, projgrid)
    
    return glats, glons, cgseas, cgpre, pretot, cgnmons, cgtcov, years, mons, cgtmp, cgprjs

def makeStacks(area, glons, glats, cgnmons, cgtcov, cgseas):    
    #this is only set up to make stacks for winter right now.
    x, y = np.meshgrid(glons, glats)
    spatialmask = ((y >= area[0,0]) & (y < area[0,1]))&((x >= area[1,0]) & (x < area[1,1]))
    
    #this is indexed to winter (winter is 0)
    cgnmons_sub = cgnmons[:, :, 0][spatialmask]
    cgtcov_sub = cgtcov[:, :, 0][spatialmask]
    cgseas_sub = cgseas[:, :, :, 0][spatialmask]
    #these thresholds can be changed. right now this just means blend all datasets
    #with more than 5 years of data, no time coverage requirement.
    inds = np.where((cgnmons_sub >= 5) & (cgtcov_sub >= 0))[0]
    
    #get anomalies
    cgseas_sub = np.divide(np.subtract(cgseas_sub.T, np.nanmean(cgseas_sub, axis = 1)), 
                           np.nanstd(cgseas_sub, axis = 1)).T
    stack = np.nanmean(cgseas_sub[inds, :], axis = 0)
    stack_std = np.nanstd(cgseas_sub[inds, :], axis = 0)
    
    return stack, stack_std

#for calculating characteristic length scales
def distance(Lats, Lons, lat, lon):
    d = np.zeros(len(Lats))
    #convert to radians    
    latsetrad = np.radians(Lats)
    lonsetrad = np.radians(Lons)
    lat = np.radians(lat)
    lon = np.radians(lon)
    
    dlat = lat- latsetrad
    dlon = lon - lonsetrad
    a = np.square(np.sin(dlat/2.0)) + np.cos(lat) * np.cos(latsetrad) * np.square(np.sin(dlon/2.0))
    #not really distance, central angle 
    d[:] = 2 * np.arcsin(np.sqrt(a))
    
    return d

def fitSemiVariogram(val, lat, lon, function):
    
    #introduce functions:
    # p[0] = c, p[1] = a
    def nugget(h, c, a):
        return np.where(h == 0, 0, c)
    def spherical(h, c, a):
        return np.where(h <= a, (c*(np.subtract(1.5*(h/a), 0.5*np.power(h/a, 3)))), c)
    def exponential(h, c, a):
        return np.where(h > 0, c*(1-np.exp(-3*h/a)), 0)
    def gaussian(h, c, a):
        return np.where(h > 0, c*(1-np.exp(-3*(np.square(h))/(a**2))), 0)
    def power(h, c, w):
        return c*np.power(h, w)
    
    if function == 'nugget':
        func = nugget
    elif function == 'spherical':
        func = spherical
    elif function == 'exponential':
        func = exponential
    elif function == 'gaussian':
        func = gaussian
    elif function == 'power':
        func = power
    else:
        print('function not recognized')
    
    maxdist = 1.5
    d = np.zeros((len(lat), len(lon)))
    davg = np.zeros(len(lat))
    disoavg = copy.copy(davg)
    diso = copy.copy(d)
    for i in range(len(lat)):
        #this is the great circle distance in radians
        d[i, :] = distance(lat, lon, lat[i], lon[i])
        #include all sites closer than the set maximum distance
        davg[i] = np.average(d[i, :][d[i, :]<= maxdist])
        #difference in isotope value
        diso[i, :] = val - val[i]
        #include all sites closer than the set maximum distance
        disoavg[i] = np.average(diso[i, :][d[i, :]<= maxdist])
    
    #semivariogramx
    svgx = []
    #semivariogramy
    svgy = []
    for i in range(len(lat)):
        svgx = np.hstack((svgx, d[i, i:]))
        svgy = np.hstack((svgy, diso[i, i:]))
    
    #establish bins 
    #what is furthest distance we want to calculate to? 1.5 rads? 
    binnum = 30
    binedges = np.arange(0, np.pi+(np.pi/binnum), np.pi/binnum)
    bins = np.zeros(binnum)
    binvar = np.zeros(binnum)
    binn = np.zeros(binnum)
    for i in range(len(bins)):
        bininds = np.where((svgx>= binedges[i])&(svgx < binedges[i+1]))[0]

        binn[i] = len(bininds)
        #variance of a finite population
        binvar[i] = 0.5*(np.nansum(np.square(svgy[bininds]))/binn[i])
        bins[i] = (binedges[i]+binedges[i+1])/2
 #-------------------------------------------------------------------------   

    nugget = 0
    traininds = np.where((~np.isnan(binvar)) & (bins < 2))[0]  
    try:
        #do not do a weighted curve fit
        popt, pcov = sp.optimize.curve_fit(func, bins[traininds], binvar[traininds], 
              bounds = (np.array([0, 0]), 
             np.array([ max(binvar), np.pi,])), 
             method = 'trf')
        pvar = np.sqrt(np.diag(pcov))
    except:
        print('problem with curve fit')
        
    #create weights matrix from semivariogram and distance matrix
    #d is the matrix with all of the distances of points from one another...
    #

    return popt, pvar, davg, disoavg

def getSemivariograms(regresults, MODregresults, MODSUBregresults, SurfEst, 
                      MODSurfEst, MODSUBSurfEst, glats, glons):
    columns = ['DJF SH 5x5 agg data', 'DJF NH 5x5 agg data', 'DJF 5x5 agg data',
               'DJF SH 5x5 ens', 'DJF NH 5x5 ens', 'DJF 5x5 ens',
               'DJF NH 5x5 ens sub', 'DJF SH 5x5 ens sub', 'DJF 5x5 ens sub', 
               'JJA NH 5x5 agg data', 'JJA SH 5x5 agg data', 'JJA 5x5 agg data',
               'JJA SH 5x5 ens', 'JJA NH 5x5 ens', 'JJA 5x5 ens',
               'JJA SH 5x5 ens sub', 'JJA NH 5x5 ens sub', 'JJA 5x5 ens sub']
    semidata = pd.DataFrame(columns = columns)
    timename = {'JJA NH':(1, 'pos'), 'JJA SH':(1, 'neg'), 'DJF NH':(0, 'pos'), 
                'DJF SH':(0, 'neg'), 'JJA':(1, 'both'), 'DJF':(0, 'both')}
    
    latinds = range(regresults.shape[0])
    loninds = range(regresults.shape[1])
    indx, indy = np.meshgrid(loninds, latinds)
    
    for i in timename.keys():
        mask = ~np.isnan(regresults[:, :, 2, timename[i][0]])
        diff = regresults[:, :, 2, timename[i][0]][mask]
        avg = SurfEst[:, :, timename[i][0]][mask]
        timeval = {'avg':avg, 'diff':diff}
        ix = glons[indx[mask]]
        iy = glats[indy[mask]]
        
        if timename[i][1] == 'pos':
            lm = iy >=0
        elif timename[i][1] == 'neg':
            lm = iy <=0  
        else:
            #this should get them all 
            lm = iy > -900
            
        popt_array = []
        for tv in timeval.keys():
            popt, pvar, davg, disoavg = fitSemiVariogram(timeval[tv][lm], iy[lm], ix[lm], 'spherical')
            popt_array.append(popt[1])
        semidata['{} 5x5 agg data'.format(i)] = np.array(popt_array)
      
        mask = ~np.isnan(MODregresults[:, :, 2, timename[i][0]])
        diff = MODregresults[:, :, 2, timename[i][0]][mask]
        avg = MODSurfEst[:, :, timename[i][0]][mask]
        timeval = {'avg':avg, 'diff':diff}
        ix = glons[indx[mask]]
        iy = glats[indy[mask]]
        
        if timename[i][1] == 'pos':
            lm = iy >=0
        elif timename[i][1] == 'neg':
            lm = iy <=0  
        else:
            #this should get them all 
            lm = iy > -900
        
        popt_array = []
        for tv in timeval.keys():
            popt, pvar, davg, disoavg = fitSemiVariogram(timeval[tv][lm], iy[lm], ix[lm], 'spherical')
            popt_array.append(popt[1])
        semidata['{} 5x5 ens'.format(i)] = np.array(popt_array)
        
        mask = ~np.isnan(MODSUBregresults[:, :, 2, timename[i][0]])
        diff = MODSUBregresults[:, :, 2, timename[i][0]][mask]
        avg = MODSUBSurfEst[:, :, timename[i][0]][mask]
        timeval = {'avg':avg, 'diff':diff}
        ix = glons[indx[mask]]
        iy = glats[indy[mask]]
        
        if timename[i][1] == 'pos':
            lm = iy >=0
        elif timename[i][1] == 'neg':
            lm = iy <=0  
        else:
            #this should get them all 
            lm = iy > -900
            
        popt_array = []
        for tv in timeval.keys():
            popt, pvar, davg, disoavg = fitSemiVariogram(timeval[tv][lm], iy[lm], ix[lm], 'spherical')
            popt_array.append(popt[1])
        semidata['{} 5x5 ens sub'.format(i)] = np.array(popt_array)
        
    return semidata

def MoransIweights(glats, glons, popt):
    def nugget(h, c, a):
        return np.where(h == 0, 0, c)
    def spherical(h, c, a):
        return np.where(h <= a, (c*(np.subtract(1.5*(h/a), 0.5*np.power(h/a, 3)))), c)
    def exponential(h, c, a):
        return np.where(h > 0, c*(1-np.exp(-3*h/a)), 0)
    def gaussian(h, c, a):
        return np.where(h > 0, c*(1-np.exp(-3*(np.square(h))/(a**2))), 0)
    def power(h, c, w):
        return c*np.power(h, w)
   
    
    lats, lons = np.meshgrid(glats, glons)    
    d = distance(lats.flatten(), lons.flatten(), lats[0,0], lons[0,0])
    vals = spherical(d, popt[0], popt[1])
    wvals = (vals.max()-vals)/vals.max()
    return wvals

def MoransI(regresults, glons, glats):
    def spherical(h, c, a):
        return np.where(h <= a, (c*(np.subtract(1.5*(h/a), 0.5*np.power(h/a, 3)))), c)

    latinds = range(regresults.shape[0])
    loninds = range(regresults.shape[1])
    indx, indy = np.meshgrid(loninds, latinds)
    timename = {'JJA':1, 'DJF':0}
    
    for i in timename.keys():
        #mask so you only get regression values that are non nans
        mask = ~np.isnan(regresults[:, :, 2, timename[i]])
        #apply the mask
        diff = regresults[:, :, 2, timename[i]][mask]
        #apply the mask to the latitudes and longitudes (these are grids) 
        ix = glons[indx[mask]]
        iy = glats[indy[mask]]  
        
        #d is the distnace between a grid cell and all other grid cells
        #diso is the difference between that grid cell value and all other grid cells
        #initalize these matrixes
        d = np.zeros((len(ix), len(iy)))
        
        for i in range(len(ix)):
            #this is the great circle distance in radians
            #these are the distances that are used to calculate the weights (w_ij) in the equation
            d[i, :] = distance(iy, ix, iy[i], ix[i])
            
        #difference of slope from average slope for the season
        #These are z's in the equation
        diso = diff - np.nanmean(diff)

        #convert distances to degrees from radians to avoid values less than 1
        d_deg = d*180/np.pi
        #take the inverse and square it
        d_deg = np.power(np.divide(1, d_deg), 2)
        #apply a search neighborhood that's the same as what we use for the semivariogram
        wvals = np.where(d > 0.5, 0, d_deg)
        wvals = np.where(np.isinf(wvals), 0, wvals)
        
        #I had used the semivariogram for the weighting... but I dont know that I like this
#        popt, pvar, davg, disoavg = fitSemiVariogram(diff, iy, ix, 'spherical')
#        vals = spherical(d, popt[0], popt[1])
#        #invert
#        wvals = (vals.max()-vals)/vals.max()
#        wvals = np.where(wvals == 1, 0, wvals)
#        #rescale
#        wvals = (wvals-wvals.min())/(wvals.max()-wvals.min())        
        
        sumsi = np.zeros(wvals.shape[0])  
        sumsw = np.zeros(wvals.shape[0])
        denom = np.sum(np.power(diso, 2))
        #svar = np.var(diff)
        keep = []
        for i in range(wvals.shape[0]):
            if (any(wvals[i, 1:]) > 0) and (i != wvals.shape[0]):
                #calculate the cross product of z_i with z_j for all j 
                zi = diso[i]
                zj = diso[(i+1):]
                temp =  np.multiply(zi, zj)
                #this is the weighting and summation of all j 
                w = wvals[i, (i+1):]
                sumsw[i] = np.sum(w)
                sumsi[i] = np.sum(np.multiply(temp, w))
                #some i don't have any weights that are non-zero
                keep.append(i)
                
        sumsi = sumsi[keep]  
        S0 = np.sum(sumsw)
        sums = np.sum(sumsi)*len(diso)/(S0*denom)
        print(sums)
        
        #calculate the z score
        #null hypothesis
        n = len(diff)
        h0 = 1/(len(diff)-1)
        print(h0)
        S0=np.sum(sumsw)
        #this matrix is symmetric wrt to weights
        S1_add = np.zeros(wvals.shape[0])
        S2_add = np.zeros(wvals.shape[0])
        for i in range(wvals.shape[0]):
            S1_add[i] = np.sum(np.power(np.add(wvals[i, (i+1):], wvals[(i+1):, i]), 2))
            S2_add[i] = np.power(np.add(np.sum(wvals[i, (i+1):]), np.sum(wvals[(i+1):, i])), 2)
            
        S1 = 0.5*np.sum(S1_add)
        
        S2 = np.sum(S2_add)
        A = n*(((n**2)+3*n+3)*S1-(n*S2)+3*(S0**2))
        z = diso
        D = np.sum(np.power(z, 4))/np.power(np.sum(np.power(z, 2)), 2)        
        B = D*((n**2 - n)*S1-2*n*S2+6*S0**2)
        C = (n-1)*(n-2)*(n-3)*S0**2
        EI_sq = (A-B)/C
        sqrtVI = np.sqrt(EI_sq-h0**2)
        zscore = (sums - h0)/sqrtVI
        print(zscore)
            
def GearysC(regresults, glons, glats):
    def spherical(h, c, a):
        return np.where(h <= a, (c*(np.subtract(1.5*(h/a), 0.5*np.power(h/a, 3)))), c)

    latinds = range(regresults.shape[0])
    loninds = range(regresults.shape[1])
    indx, indy = np.meshgrid(loninds, latinds)
    timename = {'JJA':1, 'DJF':0}
    
    for i in timename.keys():
        mask = ~np.isnan(regresults[:, :, 2, timename[i]])
        diff = regresults[:, :, 2, timename[i]][mask]
        ix = glons[indx[mask]]
        iy = glats[indy[mask]]       
        #d is the distnace between a grid cell and all other grid cells
        #diso is the difference between that grid cell value and all other grid cells

        d = np.zeros((len(ix), len(iy)))
        diso = copy.copy(d)
        
        for i in range(len(ix)):
            #this is the great circle distance in radians
            d[i, :] = distance(iy, ix, iy[i], ix[i])
            #difference of slope at one point from all other points
            diso[i, :] = diff[i] - diff
            
        d_deg = d*180/np.pi
        #take the inverse and square it
        d_deg = np.power(np.divide(1, d_deg), 2)
        #apply a search neighborhood that's the same as what we use for the semivariogram
        wvals = np.where(d > 1, 0, d_deg)
        
#        popt, pvar, davg, disoavg = fitSemiVariogram(diff, iy, ix, 'spherical')
#        vals = spherical(d, popt[0], popt[1])
#        #invert
#        wvals = (vals.max()-vals)/vals.max()
#        wvals = np.where(wvals == 1, 0, wvals)
#        #rescale
#        wvals = (wvals-wvals.min())/(wvals.max()-wvals.min())
        
        sumsi = np.zeros(wvals.shape[0])
        sumsw = np.zeros(wvals.shape[0])
        denom = np.zeros(wvals.shape[0])
        keep = []
        for i in range(wvals.shape[0]):
            if (any(wvals[i, 1:]) > 0) and (i != wvals.shape[0]):
                sumsi[i] =  np.sum(np.multiply(wvals[i, (i+1):], np.power(diso[i, (i+1):], 2)))
                sumsw[i] = np.sum(wvals[i, (i+1):])
                denom[i] = np.power(diff[i]-np.average(diff), 2)
                keep.append(i)
                
        sumsi = sumsi[keep]  
        sums = (len(diff)-1)*np.sum(sumsi)/(2*np.sum(sumsw)*np.sum(denom))
        print(sums)        
        
        
        
        
        


        
    
    
    