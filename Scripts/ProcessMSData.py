# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:39:30 2019

@author: u0929173
"""

def processCAM(os, np, nc4, copy, GetData, CAMdir):
    filename = """f.e13.FiPIPDC5CN.5x5.ctrl.cam.h0.194001-201412.nc"""
    (Lons, Lats, icmons, 
        icyears, dates, d18O0, 
        Pre0, Tmp0, PW0, grids0, 
        pregrids0, tmpgrids0, pwgrids0) = GetData.getCAM(os, np, nc4, copy, filename, CAMdir)
        
    filename = """f.e13.FiPIPDC5CN.5x5.ens1.cam.h0.194001-201412.nc"""
    (Lons, Lats, icmons, 
        icyears, dates, d18O1, 
        Pre1, Tmp1, PW1, grids1, 
        pregrids1, tmpgrids1, pwgrids1) = GetData.getCAM(os, np, nc4, copy, filename, CAMdir)
    
    filename = """f.e13.FiPIPDC5CN.5x5.ens2.cam.h0.194001-201412.nc"""
    (Lons, Lats, icmons, 
        icyears, dates, d18O2, 
        Pre2, Tmp2, PW2, grids2, 
        pregrids2, tmpgrids2, pwgrids2) = GetData.getCAM(os, np, nc4, copy, filename, CAMdir)
    
    d18O =  np.nanmean(np.concatenate((d18O0[None], d18O1[None], d18O2[None]), axis = 0), axis = 0)   
    Pre =  np.nanmean(np.concatenate((Pre0[None], Pre1[None], Pre2[None]), axis = 0), axis = 0) 
    Tmp =  np.nanmean(np.concatenate((Tmp0[None], Tmp1[None], Tmp2[None]), axis = 0), axis = 0) 
    PW =  np.nanmean(np.concatenate((PW0[None], PW1[None], PW2[None]), axis = 0), axis = 0) 
    
    grids = np.nanmean(np.concatenate((grids0[None], grids1[None], grids2[None]), axis = 0), axis = 0)
    pregrids = np.nanmean(np.concatenate((pregrids0[None], pregrids1[None], pregrids2[None]), axis = 0), axis = 0)  
    tmpgrids = np.nanmean(np.concatenate((tmpgrids0[None], tmpgrids1[None], tmpgrids2[None]), axis = 0), axis = 0)  
    pwgrids = np.nanmean(np.concatenate((pwgrids0[None], pwgrids1[None], pwgrids2[None]), axis = 0), axis = 0) 
    
    return Lons, Lats, icyears, d18O, Pre, Tmp, PW, grids, pregrids, tmpgrids, pwgrids

def processData(os, np, nc4, csv, copy, GetData):
    (t, glons, glats, coarsegrid, 
     pregrid, gridDelT, ninds) = GetData.getdata(os, np, nc4, csv, copy)

    (cgseas, cgpre, pretot, 
     cgnmons, cgtcov, years, mons) = GetData.getseasonalgrids(np, glats, glons, 
                                        coarsegrid, pregrid, t)
    
    return glats, glons, cgseas, cgpre, pretot, cgnmons, cgtcov, years, mons

def makeStacks(np, area, glons, glats, cgnmons, cgtcov, cgseas):    
    #this is only set up to make stacks for winter right now.
    x, y = np.meshgrid(glons, glats)
    spatialmask = ((y >= area[0,0]) & (y < area[0,1]))&((x >= area[1,0]) & (x < area[1,1]))
    
    #this is indexed to winter (summer is 0)
    cgnmons_sub = cgnmons[:, :, 1][spatialmask]
    cgtcov_sub = cgtcov[:, :, 1][spatialmask]
    cgseas_sub = cgseas[:, :, :, 1][spatialmask]
    #these thresholds can be changed. right now this just means blend all datasets
    #with more than 5 years of data, no time coverage requirement.
    inds = np.where((cgnmons_sub >= 5) & (cgtcov_sub >= 0))[0]
    
    #get anomalies
    cgseas_sub = np.divide(np.subtract(cgseas_sub.T, np.nanmean(cgseas_sub, axis = 1)), np.nanstd(cgseas_sub, axis = 1)).T
    stack = np.nanmean(cgseas_sub[inds, :], axis = 0)
    stack_std = np.nanstd(cgseas_sub[inds, :], axis = 0)
    
    return stack, stack_std

#for calculating characteristic length scales
def distance(np, Lats, Lons, lat, lon):
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

def fitSemiVariogram(os, sp, np, mpl, copy, distance, val, lat, lon, 
                     function):
    
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
        d[i, :] = distance(np, lat, lon, lat[i], lon[i])
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
        print 'problem with curve fit'
#    viridis = mpl.pyplot.get_cmap('viridis')
#    mpl.pyplot.figure()
#    mpl.pyplot.scatter(bins, binvar, 30, c = binn, cmap = viridis, linewidth = 0.5)
#    mpl.pyplot.scatter(bins[traininds], binvar[traininds], c = binn[traininds], cmap = viridis, edgecolors = 'k', linewidth = 0.5)
#    mpl.pyplot.colorbar()
#    mpl.pyplot.plot(bins, func(bins, popt[0], popt[1]), color = 'k', linewidth = 2)
#    mpl.pyplot.plot([popt[1], popt[1]], [0, max(binvar)*1.10], color = 'gray', linewidth = 1.5)   
#    mpl.pyplot.plot([0, 3], [popt[0], popt[0]], color = 'gray', linewidth = 1.5)                 
#    mpl.pyplot.xlabel('Lag distance, radians')
#    mpl.pyplot.ylabel('Semivariance')
    
    return popt, pvar, davg, disoavg

def getSemivariograms(os, np, mpl, sp, pd, copy, distance, fitSemiVariogram, regresults, 
                      MODregresults, MODSUBregresults, SurfEst, MODSurfEst, 
                      MODSUBSurfEst, glats, glons):
    columns = ['DJF 5x5 agg data', 'DJF 5x5 ens',  'DJF 5x5 ens sub', 
               'JJA 5x5 agg data', 'JJA 5x5 ens', 'JJA 5x5 ens sub']
    semidata = pd.DataFrame(columns = columns)
    timename = {'JJA':0, 'DJF':1}
    
    latinds = range(regresults.shape[0])
    loninds = range(regresults.shape[1])
    indx, indy = np.meshgrid(loninds, latinds)
    
    for i in timename.keys():
        mask = ~np.isnan(regresults[:, :, 2, timename[i]])
        diff = regresults[:, :, 2, timename[i]][mask]
        avg = SurfEst[:, :, timename[i]][mask]
        timeval = {'avg':avg, 'diff':diff}
        ix = glons[indx[mask]]
        iy = glats[indy[mask]]
        
        popt_array = []
        for tv in timeval.keys():
            popt, pvar, davg, disoavg = fitSemiVariogram(os, sp,
                                  np, mpl, copy, distance, timeval[tv], iy, ix, 'spherical')
            popt_array.append(popt[1])
        semidata['{} 5x5 agg data'.format(i)] = np.array(popt_array)
      
        mask = ~np.isnan(MODregresults[:, :, 2, timename[i]])
        diff = MODregresults[:, :, 2, timename[i]][mask]
        avg = MODSurfEst[:, :, timename[i]][mask]
        timeval = {'avg':avg, 'diff':diff}
        ix = glons[indx[mask]]
        iy = glats[indy[mask]]
        popt_array = []
        for tv in timeval.keys():
            popt, pvar, davg, disoavg = fitSemiVariogram(os, sp,
                                  np, mpl, copy, distance, timeval[tv], iy, ix, 'spherical')
            popt_array.append(popt[1])
        semidata['{} 5x5 ens'.format(i)] = np.array(popt_array)
        
        mask = ~np.isnan(MODSUBregresults[:, :, 2, timename[i]])
        diff = MODSUBregresults[:, :, 2, timename[i]][mask]
        avg = MODSUBSurfEst[:, :, timename[i]][mask]
        timeval = {'avg':avg, 'diff':diff}
        ix = glons[indx[mask]]
        iy = glats[indy[mask]]
        popt_array = []
        for tv in timeval.keys():
            popt, pvar, davg, disoavg = fitSemiVariogram(os, sp,
                                  np, mpl, copy, distance, timeval[tv], iy, ix, 'spherical')
            popt_array.append(popt[1])
        semidata['{} 5x5 ens sub'.format(i)] = np.array(popt_array)
        
    return semidata