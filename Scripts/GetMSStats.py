# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:13:05 2019

@author: u0929173
"""
import numpy as np
import statsmodels.api as sm
import scipy.odr as odr
import matplotlib as mpl

def areaPerGridCell(glats, glons):
    #radius of earth
    latareas = np.zeros(len(glats))
    r=6378.1
    lonboxes = len(glons)
    lat0=np.pi/2
    for i, lat in enumerate(np.flip(glats[18:])):        
        #latitude in radians
        lat = (np.pi/180)*lat  
        print(lat0, lat)
        A1 = (2*np.pi*r**2)*(np.sin(lat0) - np.sin(lat))/lonboxes
        latareas[i] = A1
        lat0=lat
    #do somthing special for the equator! since bounding boxes are +2.5 and -2.5
    A1=2*(2*np.pi*r**2)*(np.sin(lat0) - np.sin(0))/lonboxes
    i=i+1
    latareas[i]=A1
    #bc symmetric     
    latareas[(i+1):] = np.flip(latareas[:(i-1)])
    totalarea = sum(latareas*lonboxes)
    latareas = np.repeat(latareas[:, None], lonboxes, axis = 1)
    return latareas, totalarea
    

def linreg(B, x):
    #local meteoric water line function
    #for use in ODR regression
    return B[0]*x +B[1]

def MaritimeSummary(df):
    maritime = df[(df['Maritime'] == 1)&(abs(df['Lat'])<= 15)]
    subset = df[(df['Maritime'] == 1)&(abs(df['Lat'])> 35)]
    print('Results for d18O slopes vs Pressure slopes:')
    y = subset['Slope'].values.astype(float)
    x = 10*subset['Press'].values.astype(float)
    x = sm.add_constant(x)
    results = sm.OLS(y, x, missing = 'drop').fit()
    print(results.summary())
      
    print('Results for d18O slopes vs Omega 500 slopes:')
    y = maritime['Slope'].values.astype(float)
    x = 10*maritime['Omega'].values.astype(float)
    x = sm.add_constant(x)
    results = sm.OLS(y, x, missing = 'drop').fit()
    print(results.summary())
    
def summary(RTregresults, regresults, latareas, totalarea):
    #summer
    seas = ['JJA', 'DJF']
    for i in range(2):
        print('----------------------------------------------')
        print('Season: {}'.format(seas[i]))
        RTslope = RTregresults[:, :, 2, i]
        DATslope = regresults[:, :, 2, i]
        mask = ~np.isnan(RTslope)&~np.isnan(DATslope)
        areas = latareas[mask]
        ss = (((RTslope[mask] > 0)*(DATslope[mask] > 0))+((RTslope[mask] < 0)*(DATslope[mask] < 0))).astype(int)
        A = sum(np.multiply(ss,areas))
        posslope = (DATslope[mask] >0).astype(int)
        posslopeA = sum(np.multiply(posslope, areas))
        print('Total # grid cells: {}\nSame sign # grid cells: {}'.format(len(ss),sum(ss))) 
        print('Total area: {}\narea with same sign: {}'.format(sum(areas), A)) 
        print('Positive Slope Area: {}'.format(posslopeA))
        perctrue = sum(ss)/float(len(ss))
        
        print('% same sign: {}'.format(perctrue*100))
        print('% of possible area: {}'.format(100*(A/sum(areas))))
        print('% of possible area is positive: {}'.format(posslopeA/sum(areas)))
        print('----------------------------------------------')
        
def comparesignplot(RTregresults, regresults):
    seas = ['JJA', 'DJF']
    for i in range(2):
        print('----------------------------------------------')
        print('Season: {}'.format(seas[i]))
        RTslope = np.where(~np.isnan(RTregresults[:, :, 2, i])&(RTregresults[:, :, 2, i] > 0), 1, RTregresults[:, :, 2, i])
        RTslope = np.where(~np.isnan(RTslope)& (RTslope < 0), -1, RTslope)
        DATslope = np.where(~np.isnan(regresults[:, :, 2, i])&(regresults[:, :, 2, i] > 0), 1, regresults[:, :, 2, i])
        DATslope = np.where(~np.isnan(DATslope)& (DATslope < 0), -1, DATslope)
        mpl.pyplot.figure()
        mpl.pyplot.pcolor(np.add(RTslope, DATslope))
        mpl.pyplot.colorbar()
        
def regNums(regresults):
    seas = ['JJA', 'DJF']
    for i in range(2):
        print('----------------------------------------------')
        print('Season: {}'.format(seas[i]))        
        #get number of regressions in the season
        nregs = np.sum(~np.isnan(regresults[:, :, 2, i]))
        nsigregs = np.sum((~np.isnan(regresults[:, :, 4, i]))&(regresults[:, :, 4, i] < 0.1))
        print('Total # regressions: {}'.format(nregs))
        print('# sig regressions (p < 0.1): {}'.format(nsigregs))
        print('----------------------------------------------')
    print('----------------------------------------------')
    print('Both') 
    
    nregs = np.sum(~np.isnan(regresults[:, :, 2, 0])&~np.isnan(regresults[:, :, 2, 1]))
    nsigregs = np.sum((~np.isnan(regresults[:, :, 4, 0]))&(regresults[:, :, 4, 0] < 0.1)&
                (~np.isnan(regresults[:, :, 4, 1]))&(regresults[:, :, 4, 1] < 0.1))
    print('Total # regression grid cells in both seasons: {}'.format(nregs))
    print('# sig regressions (p < 0.1) in both seasons: {}'.format(nsigregs))
    print('----------------------------------------------')
    
def dataNums(regresults, cgnmons):
    seas = ['JJA', 'DJF']
    for i in range(2):
        print('----------------------------------------------')
        print('Season: {}'.format(seas[i]))        
        #get number of regressions in the season
        nregs = np.sum(cgnmons[:, :, i][~np.isnan(regresults[:, :, 2, i])])
        print('Total # grid-seasons of data: {}'.format(nregs))
        print('----------------------------------------------')    

def regressSlopes(regresults, MODregresults):
    print('Results for d18O slopes vs Pressure slopes:')
    mask = ~np.isnan(regresults[:, :, 2, :])
    y = regresults[:, :, 2, :][mask]
    x = MODregresults[:, :, 2, :][mask]
    corr = np.corrcoef(x, y)
    print('Correlation coefficient: {}'.format(corr[0,1]))
    x = sm.add_constant(x)
    results = sm.OLS(y, x).fit()
    print(results.summary())

def MaritimeSummary_ODR(df):
    
    maritime = df[(df['Maritime'] == 1)&(abs(df['Lat'])<= 15)]
    subset = df[(df['Maritime'] == 1)&(abs(df['Lat'])> 35)]
    
    linmod = odr.Model(linreg)
    #instantiate with results of linear regression
    beta0 = [0.0012, 0.0059]
    y = subset['Slope'].values.astype(float)
    x = 10*subset['Press'].values.astype(float)
    sx = subset['Slope_std'].values.astype(float)
    sy = subset['Press_std'].values.astype(float)
    realdata = odr.RealData(x, y, sx=sx, sy=sy)
    model = odr.ODR(realdata, linmod, beta0 = beta0)
    results = model.run()
    print('Beta0: {}, Beta0 std: {}, Beta1: {}, Beta1 std: {}, inv condnum: {}, resvar: {}'.format(results.beta[0], results.sd_beta[0], 
                                         results.beta[1], results.sd_beta[1], 
                                         results.inv_condnum, results.res_var))

    x = sm.add_constant(x)
    results = sm.OLS(y, x, missing = 'drop').fit()
    print(results.summary())
      
    print('Results for d18O slopes vs Omega 500 slopes:')
    y = maritime['Slope'].values.astype(float)
    x = 10*maritime['Omega'].values.astype(float)
    x = sm.add_constant(x)
    results = sm.OLS(y, x, missing = 'drop').fit()
    print(results.summary())