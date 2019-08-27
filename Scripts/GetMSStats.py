# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:13:05 2019

@author: u0929173
"""

def MaritimeSummary(np, sm, df):
    maritime = df[(df['Maritime'] == 1)&(abs(df['Lat'])<= 20)]
    subset = df[(df['Maritime'] == 1)&(abs(df['Lat'])> 25)]
    print('Results for d18O slopes vs Pressure slopes:')
    y = subset['Slope'].values
    x = 10*subset['Press'].values.astype(float)
    x = sm.add_constant(x)
    results = sm.OLS(y, x).fit()
    print(results.summary())
      
    print('Results for d18O slopes vs Omega 500 slopes:')
    y = maritime['Slope'].values
    x = 10*maritime['Omega'].values.astype(float)
    x = sm.add_constant(x)
    results = sm.OLS(y, x).fit()
    print(results.summary())
    
def summary(np, RTregresults, regresults):
    #summer
    seas = ['JJA', 'DJF']
    for i in range(2):
        print('----------------------------------------------')
        print('Season: {}'.format(seas[i]))
        RTslope = RTregresults[:, :, 2, i]
        DATslope = regresults[:, :, 2, i]
        mask = ~np.isnan(RTslope)&~np.isnan(DATslope)
        ss = (((RTslope[mask] > 0)*(DATslope[mask] > 0))+((RTslope[mask] < 0)*(DATslope[mask] < 0))).astype(int)
        print('Total # grid cells: {}\nSame sign # grid cells: {}'.format(len(ss),sum(ss)))
        perctrue = sum(ss)/float(len(ss))
        print('% same sign: {}'.format(perctrue*100))
        print('----------------------------------------------')

def regNums(np, regresults):
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
    
def dataNums(np, regresults, cgnmons):
    seas = ['JJA', 'DJF']
    for i in range(2):
        print('----------------------------------------------')
        print('Season: {}'.format(seas[i]))        
        #get number of regressions in the season
        nregs = np.sum(cgnmons[:, :, i][~np.isnan(regresults[:, :, 2, i])])
        print('Total # grid-seasons of data: {}'.format(nregs))
        print('----------------------------------------------')    

def regressSlopes(np, sm, regresults, MODregresults):
    print('Results for d18O slopes vs Pressure slopes:')
    mask = ~np.isnan(regresults[:, :, 2, :])
    y = regresults[:, :, 2, :][mask]
    x = MODregresults[:, :, 2, :][mask]
    corr = np.corrcoef(x, y)
    print('Correlation coefficient: {}'.format(corr[0,1]))
    x = sm.add_constant(x)
    results = sm.OLS(y, x).fit()
    print(results.summary())
    