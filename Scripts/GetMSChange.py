# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:09:33 2019

@author: u0929173
"""

def GetMultidecadalChange(np, sm, glats, glons, years, cgseas, cgtcov, cgnmons):   
    regresults = np.zeros((len(glats), len(glons), 6, 2))
    regresults.fill(np.nan)
    
    for i in range(cgseas.shape[0]):
        for j in range(cgseas.shape[1]):
            for k in range(cgseas.shape[3]):
                if (cgnmons[i, j, k]>= 10) & (cgtcov[i, j, k] >= 25):
                    y = cgseas[i, j, :, k]
                    mask = ~np.isnan(y)
                    y = y[mask]                    
                    #make it so that it's calculating change per decade!
                    x = np.unique(years)[mask]/10
                    x = sm.add_constant(x)
                    # if you want to analyze the anomaly
                    #y = (y-np.average(y))/np.std(y)
                    #x = sm.add_constant((x - np.average(x))/np.std(x))                    
                    results = sm.OLS(y, x).fit()
                    
                    savevals = np.array([results.params[0], results.bse[0], results.params[1], results.bse[1], results.f_pvalue, results.rsquared])
                    regresults[i, j, :, k] = savevals
                    
                    
    return regresults

def getCAMChange(np, sm, icyears, years, Lats, Lons, grids, cgnmons, cgtcov, thresh):
    #get the regression based on same times/locations as data. Use data as a mask.
    startind = np.where(np.unique(icyears) == 1960)[0][0]
    yearend = np.where(np.unique(years) == np.max(np.unique(icyears[startind:])))[0][0]+1   
    icsregresults = np.zeros((len(Lats), len(Lons), 6, 2))
    icsregresults.fill(np.nan)
    
    for i in range(grids.shape[0]):
        for j in range(grids.shape[1]):
            for k in range(grids.shape[3]):
                if (cgnmons[i, j, k]>= thresh[0])& (cgtcov[i, j, k] >= thresh[1]):
                    y = grids[i, j, startind:, k]
                    mask = ~np.isnan(y)
                    y = y[mask]
                    if len(y) > 5:
                        x = sm.add_constant((np.unique(icyears)/10)[startind:][mask])
                        results = sm.OLS(y, x).fit()
                        savevals = np.array([results.params[0], results.bse[0], results.params[1], results.bse[1], results.f_pvalue, results.rsquared])
                        icsregresults[i, j, :, k] = savevals 
        
                    else:
                        savevals = np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
                        icsregresults[i, j, :, k] = savevals  

    return icsregresults

