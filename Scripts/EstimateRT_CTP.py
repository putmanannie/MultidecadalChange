# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 10:28:11 2019

@author: u0929173
"""

def EstModRT(np, copy, icyears, grids, pregrids, pwgrids, tmpgrids):
    startyr = np.array([1960, 1960, 1995])
    endyr = np.array([2015, 1980, 2015])
    
    prediff = np.zeros((grids.shape[0], grids.shape[1], 3, 2))
    tmpdiff = copy.copy(prediff)
    pwdiff = copy.copy(prediff)
    griddiff = copy.copy(prediff)
    d18O_est = np.zeros((grids.shape[0], grids.shape[1], 2, 2))
    
    for i in range(3):
        sind = np.argmin(abs(np.unique(icyears)-startyr[i]))
        eind = np.argmin(abs(np.unique(icyears)-endyr[i]))
        prediff[:, :, i, :] = np.nanmean(pregrids[:, :, sind:eind, :], axis = 2) 
        tmpdiff[:, :, i, :] = np.nanmean(tmpgrids[:, :, sind:eind, :], axis = 2) 
        pwdiff[:, :, i, :] = np.nanmean(pwgrids[:, :, sind:eind, :], axis = 2)
        griddiff[:, :, i, :]= np.nanmean(grids[:, :, sind:eind, :], axis = 2)
        
        if i > 0:
            Qanom = np.divide(np.subtract(pwdiff[:, :, i, :], pwdiff[:, :, 0, :]), pwdiff[:, :, 0, :])
            Panom = np.divide(np.subtract(prediff[:, :, i, :], prediff[:, :, 0, :]), prediff[:, :, 0, :])
            d18O_est[:, :, (i-1), :] = np.add(np.multiply(-0.4*np.subtract(Qanom, Panom), griddiff[:, :, 0, :]), griddiff[:, :, 0, :])
    return prediff, tmpdiff, pwdiff, griddiff, d18O_est

def PlotModRT(os, np, mpl, ctp, ccrs, aggplots, Lons, Lats, d18O_est, griddiff, prediff, tmpdiff, plottingmask):
    cmap = mpl.pyplot.get_cmap('RdBu_r')  
    BrBG = mpl.pyplot.get_cmap('BrBG_r')
    
    x, y = np.meshgrid(Lons, Lats)
    d18Odiff_est = np.subtract(d18O_est[:, :, 1, :], d18O_est[:, :, 0, :])
    d18Odiff_mod = np.subtract(griddiff[:, :, 2, :], griddiff[:, :, 1, :])
    prediff_mod = np.subtract(prediff[:, :, 2, :], prediff[:, :, 1, :])
    DelT = np.subtract(tmpdiff[:, :, 2, :], tmpdiff[:, :, 1, :])
    #get all locations where both estiamtes are greater than 0
    sspos = (d18Odiff_est > 0)*(d18Odiff_mod > 0)
    #get all locations where both estimates are less than 0
    ssneg = (d18Odiff_est < 0)*(d18Odiff_mod < 0)
    
    ss = (sspos+ssneg).astype(int)
    ss = np.ma.masked_less(ss, 0.5)
    seas = ['JJA', 'DJF']
    for i in range(2):        
        mpl.pyplot.figure()
        mpl.pyplot.scatter(DelT[:, :, i][plottingmask[:, :, i]], d18Odiff_mod[:, :, i][plottingmask[:, :, i]], 
                           c = prediff_mod[:, :, i][plottingmask[:, :, i]], vmin = -1.5, vmax = 1.5, cmap = BrGB, edgecolors = 'k')
        mpl.pyplot.plot([-0.5, 4], [0, 0], c = 'gray', linewidth = 2, zorder = 0)
        mpl.pyplot.plot([0, 0], [-0.5, 1], c = 'gray', linewidth = 2, zorder = 0)
        mpl.pyplot.colorbar()
        mpl.pyplot.xlabel(" $\Delta T$, iCAM")
        mpl.pyplot.ylabel("$\Delta \delta^{18}O$, iCAM")
        mpl.pyplot.savefig(os.path.join(aggplots, 'T_P_Deld18O_{}.png'.format(seas[i])), dpi = 500)
        
        #figure plto only where we have data (need to run other script to make this work)
      
        mpl.pyplot.figure()
        p = mpl.pyplot.scatter(d18Odiff_mod[:, :, i][plottingmask[:, :, i]], 
                           d18Odiff_est[:, :, i][plottingmask[:, :, i]], 
                            c = DelT[:, :, i][plottingmask[:, :, i]], cmap = cmap, 
                            vmin = -2, vmax = 2, edgecolors = 'k')
        mpl.pyplot.plot([-1, 1], [-1, 1], c = 'gray', linewidth = 2)
        mpl.pyplot.xlabel('iCAM $\Delta \delta^{18}O$')
        mpl.pyplot.ylabel('Residence Time $\Delta \delta^{18}O$')
        cb = mpl.pyplot.colorbar(p)
        #cb = mpl.pyplot.colorbar(p, orientation = 'horizontal', cax = ax)
        cb.set_label('$\Delta T$')
        mpl.pyplot.savefig(os.path.join(aggplots, 'iCAM_ResTime_Compare_{}.png'.format(seas[i])), dpi = 500)

    #stipple where the sign is the same
    #get same sign locations
    fig = mpl.pyplot.figure(figsize = (9, 7))
    ax1 = mpl.pyplot.subplot(2, 1, 1, projection=ccrs.Robinson(central_longitude = 0, globe = None))
    ax2 = mpl.pyplot.subplot(2, 1, 2, projection=ccrs.Robinson(central_longitude = 0, globe = None))
    axes = np.array([ax1, ax2])
    x, y = np.meshgrid(Lons, Lats)
    for (ax, sn, k) in zip(axes, seas, range(len(axes))):
        #make a plot of all grids and where there is data
        ax.add_feature(ctp.feature.COASTLINE, zorder=10)

 
        vals = np.ma.array(d18Odiff_mod[:, :, k]/3.5, mask = np.isnan(d18Odiff_mod[:, :, k]))
        p = ax.pcolormesh(x, y, vals, cmap = BrBG, vmin = -0.5, vmax = 0.5, 
                      edgecolors = 'k', linewidths = 0.5, transform = ccrs.PlateCarree())  
        ax.pcolor(x, y, ss[:, :, k], hatch = '....', zorder = 5, alpha = 0.0, transform = ccrs.PlateCarree())

    axes[0].annotate('(a) JJA', xy = (0.02, 0.95), xytext = (0.02, 0.95), xycoords = 'axes fraction')    
    axes[1].annotate('(b) DJF', xy = (0.02, 0.95), xytext = (0.02, 0.95), xycoords = 'axes fraction')
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    cb = fig.colorbar(p, cax=cbar_ax)
    #cb = mpl.pyplot.colorbar(p, orientation = 'horizontal', cax = ax)
    cb.set_label('$\Delta \delta^{18}O$ ('+u'\u2030'+' decade $^{-1}$)')
    
    mpl.pyplot.savefig(os.path.join(aggplots, 'd18O_change_iCAM_RT_Stipple.png'), dpi = 500)
    
    #get percentage where sign is the same:
    #where true 
    wheretrue = sum(ss.data.flatten())
    totalcells = ss.shape[0]*ss.shape[1]*ss.shape[2]
    perctrue = wheretrue/float(totalcells)
    print('Percent of all pixels which share a sign with RT estimate: {}'.format(perctrue*100))
        
        
        
#%%%%%%%%%%%%%%%%%%%%%%%% data RT estimate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
def GetRTData(os, np, dt, urllib, prefileindir):
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
                                   """W:\Annie\Reanalysis\{}""".format(pwfn))
                urllib.urlcleanup()
                
            if prfn not in prefileindir:
                urllib.urlretrieve("""ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface_gauss/{}""".format(prfn),
                                   """W:\Annie\Reanalysis\{}""".format(prfn))
                urllib.urlcleanup()
                
            if tpfn not in prefileindir:
                urllib.urlretrieve("""ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface/{}""".format(tpfn),
                                   """W:\Annie\Reanalysis\{}""".format(tpfn))
                urllib.urlcleanup()
                
def regridRTdata(iris, np, os):
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


def getMeasRT(os, np, nc4, copy, mons, timename, SurfEst, SurfEstSE): 
    predaysonly = False
    startyr = [1965, 1965, 2000]
    endyr = [2015, 1980, 2015]
    predaysonly = False
    
    N = np.zeros(len(startyr))
    initialize = 0 
    for i in range(len(startyr)):
        time = np.array([], dtype = 'datetime64[D]')
        #make the whole record for the timeperiod then mask it/average it at the end
        for y in range(startyr[i], endyr[i]+1):       
            
            #open the precipitation rate file
            dspr = nc4.Dataset(r"""E:\WorkComputer\WDrive\Annie\Reanalysis\Regridded\regrid_prate.sfc.gauss.{}.nc""".format(y))
            #get precipitable water
            dspw = nc4.Dataset(r"""E:\WorkComputer\WDrive\Annie\Reanalysis\Regridded\regrid_pr_wtr.eatm.{}.nc""".format(y))
            if i > 0 or initialize == 0:
                tmpr = nc4.Dataset(r'E:\WorkComputer\WDrive\Annie\Reanalysis\Regridded\regrid_air.sig995.{}.nc'.format(y))
            
            if initialize == 0:
                lat = dspr.variables['lat'][:]
                lon = dspr.variables['lon'][:]
                
                prate = np.empty((0, len(lat), len(lon)))
                prwt = copy.copy(prate)
                tmp = copy.copy(prate)
                
                P = np.empty((len(startyr), len(lat), len(lon)))
                Q = copy.copy(P)
                T = copy.copy(P)
                Pse = copy.copy(P)
                Qse = copy.copy(P)
                Tse = copy.copy(P)
                initialize = 1
                
            #this is how we handle time
            #time is same in all files
            t = dspr.variables['time'][:]/24
            t = t.astype('timedelta64[D]')
            basetime = np.datetime64('1800-01-01', 'D')
            time= np.append(time, t+basetime)
            
            #get precipitation rate and close
            #this is in kg/m**2/s (or mm/s)
            prate = np.concatenate((prate, dspr.variables['prate'][:]), axis = 0)
            #this is in kg/m**2
            prwt = np.concatenate((prwt, dspw.variables['pr_wtr'][:]), axis = 0)
            if i > 0:
                tmp = np.concatenate((tmp, tmpr.variables['air'][:]), axis = 0)
                tmpr.close()
                
            dspw.close()
            dspr.close()   
        
        
        #the the months
        monmask = np.array([ti in mons for ti in time.astype('datetime64[M]').astype(int)%12+1])
        prwt = prwt[monmask]
        prate = prate[monmask]
        if i > 0:
            tmp = tmp[monmask]
        time = time[monmask]
        #if DJF, need to remove first 1, 2 and last 12
        if timename == 'DJF':
            #remove first two months
            monmask1 = np.array((time.astype('datetime64[Y]').astype(int)+1970 == startyr[i]) & 
                ((time.astype('datetime64[M]').astype(int)%12+1 == 1) | (time.astype('datetime64[M]').astype(int)%12+1 == 2)))
            #remove last month of last year
            monmask2 = np.array((time.astype('datetime64[Y]').astype(int)+1970 == endyr[i]) & 
                (time.astype('datetime64[M]').astype(int)%12+1 == 12))
            
            time = time[~monmask1*~monmask2]
            if i > 0: 
                tmp = tmp[~monmask1*~monmask2]
            prate = prate[~monmask1*~monmask2]
            prwt = prwt[~monmask1*~monmask2]
        #calculate the P and Q for the time frame in each grid cell
        #anything less than or equal to a trace of precipitation is not considered
        if predaysonly:
            prate = np.where(prate*3600*24 <= 0.1, np.nan, prate)
        else:
            prate = np.where(prate < 0, 0, prate)
            
        prate = np.where(prate < 10, prate, np.nan)
        prwt = np.where((prwt < 1000)&(prwt > 0), prwt, np.nan)
        tmp = np.where(tmp < 1000, tmp, np.nan)       
        
        N[i] = np.shape(prwt)[0]
    #    if i > 0:
    #        #calculate the anomaly
    #        prate = np.divide(np.subtract(prate, P[np.newaxis, 0, :, :]), P[np.newaxis, 0, :, :])
        if i > 0: 
            T[i, :, :] = np.nanmean(tmp, axis = 0) 
            Tse[i, :, :] = np.nanstd(tmp, axis = 0)/np.sqrt(tmp.shape[0])    
        
        P[i, :, :] = np.nanmean(prate, axis = 0)
        Pse[i, :, :] = np.nanstd(prate, axis = 0)/np.sqrt(prate.shape[0])
        #P_hi[i, :, :] = np.nanmean()    
        
        prwt = np.where(prwt < 0, np.nan, prwt)
        #select precipitable water only for days with precipitation
        if predaysonly:
            if i > 0:
                prwt = np.where(np.isnan(prate), np.nan, prwt)
    
        Q[i, :, :] = np.nanmean(prwt, axis = 0)
        #calculate the standard error of the mean (stdev/n)
        Qse[i, :, :] = np.nanstd(prwt, axis = 0)/np.sqrt(prwt.shape[0])
    
        prate = np.empty((0, len(lat), len(lon)))
        prwt = copy.copy(prate)
        tmp = copy.copy(prate)    
###################################################
    #need to flip the values we use
    flipind = np.where(lon > 180)[0][0]
    lon = np.concatenate((lon[flipind:]-360, lon[:flipind]))
    Q = np.concatenate((Q[:, :, flipind:], Q[:, :, :flipind]), axis = 2)
    Qse = np.concatenate((Qse[:, :, flipind:], Qse[:, :, :flipind]), axis = 2) 
    P = np.concatenate((P[:, :, flipind:], P[:, :, :flipind]), axis = 2)
    Pse = np.concatenate((Pse[:, :, flipind:], Pse[:, :, :flipind]), axis = 2)
    T = np.concatenate((T[:, :, flipind:], T[:, :, :flipind]), axis = 2)
    Tse = np.concatenate((Tse[:, :, flipind:], Tse[:, :, :flipind]), axis = 2)            
    
    DRT = np.zeros((len(startyr)-1, len(lat), len(lon)))
    RT = copy.copy(DRT)
    RT_se = copy.copy(RT)
    DRT_se = copy.copy(DRT)
    Dd18O = copy.copy(DRT)
    d18O_est = copy.copy(DRT)
    d18O_est_se = copy.copy(DRT)    
    
    for i in range(len(startyr)):
        if i == 0:
            barRT = np.divide(Q[i], P[i]*3600*24)
        if i > 0:
            #these are only the anomalies, not the actual numbers.
            DQ = np.divide(np.subtract(Q[i], Q[0]), Q[0])
            DP = np.divide(np.subtract(P[i], P[0]), P[0])
            DRT[(i-1), :, :] = np.subtract(DQ[i], DP[i])
            
            RT[(i-1), :, :] = np.divide(Q[i], P[i]*3600*24)
            RT_se[(i-1), :, :] = np.multiply(RT[(i-1), :, :], np.sqrt(np.add(np.power(np.divide(Qse[i], Q[i]), 2), 
                            np.power(np.divide(Pse[i]*3600*24, P[i]*3600*24), 2))))
    
            Dd18O[(i-1), :, :] = -0.4*DRT[(i-1), :, :]
            
            #propagate the standard error of the mean also
            DRT_se[(i-1), :, :] = np.sqrt(np.add(np.power(Qse[i], 2), np.power(Pse[i], 2)))
            d18O_est[(i-1), :, :] = np.add(np.multiply(Dd18O[(i-1), :, :], SurfEst), SurfEst)
            
            temp = np.multiply(np.multiply(Dd18O[(i-1), :, :], SurfEst), 
                               np.sqrt(np.add(np.power(np.divide(DRT_se[(i-1)], DRT[(i-1)]), 2), 
                                              np.power(np.divide(SurfEstSE, SurfEst), 2))))
            d18O_est_se[(i-1), :, :] = np.sqrt(np.add(np.power(temp, 2), np.power(SurfEstSE, 2)))
           
      
    denom = np.sqrt(np.add(np.square(RT_se[0, :, :]), np.square(RT_se[1, :, :])))
    tstat = np.divide(np.subtract(RT[1], RT[0]), denom)
    #take the minimum N number to calculate degrees of freedom
    
    
    tpdf = np.random.standard_t(min(N), size = 10000)
    pval = np.zeros(tstat.shape)
    for j in range(tstat.shape[0]):
        for k in range(tstat.shape[1]):
            pval[j, k] = np.sum(tpdf[ tpdf > abs(tstat[j, k])])/float(len(tpdf))
    
    timedel = ((startyr[2]+endyr[2])/2)-((startyr[1]+endyr[1])/2)
    #save the calculated values to a netcdf file
    os.chdir(r"""E:\WorkComputer\WDrive\Annie\NCfiles""")
    f = nc4.Dataset('RT_model_d18O_pred_{}.nc'.format(timename), 'w', format = 'NETCDF4')
    ntimes = f.createDimension('times', np.shape(d18O_est)[0])
    lat = f.createDimension('lat', np.shape(T)[1])
    lon = f.createDimension('lon', np.shape(T)[2])
    single = f.createDimension('single', 1)
    
    TMP = f.createVariable('tmp', 'f', ('times', 'lat', 'lon'))
    TMP[:] = T[1:3]
    
    d18OEST = f.createVariable('d18Oest', 'f', ('times', 'lat', 'lon'))
    d18OEST[:] = d18O_est
    
    td = f.createVariable('td', 'f', ('single',))
    td[:] = timedel
    f.close()
    
    return T, d18O_est, RT

def getRT4slope(os, dt, np, nc4, copy): 
    startyr = 1960
    endyr = 2016
    yrs = range(startyr, endyr+1)
    seasdict = {'JJA': (7, 0), 'DJF': (1, 1)}    
    initialize = 0 
    #first, get timeseries of t, pw, and pre
    time = []
    #make the whole record for the timeperiod then mask it/average it at the end
    for y in yrs:      
        #open the precipitation rate file
        dspr = nc4.Dataset(r"""E:\WorkComputer\WDrive\Annie\Reanalysis\Regridded\regrid_prate.sfc.gauss.{}.nc""".format(y))
        #get precipitable water
        dspw = nc4.Dataset(r"""E:\WorkComputer\WDrive\Annie\Reanalysis\Regridded\regrid_pr_wtr.eatm.{}.nc""".format(y))
        #get temperature
        tmpr = nc4.Dataset(r'E:\WorkComputer\WDrive\Annie\Reanalysis\Regridded\regrid_air.sig995.{}.nc'.format(y))
        
        if initialize == 0:
            lat = dspr.variables['lat'][:]
            lon = dspr.variables['lon'][:]

            prate = np.empty((0, len(lat), len(lon)))
            prwt = copy.copy(prate)
            tmp = copy.copy(prate)
            initialize = 1
           
        #this is how we handle time
        #time is same in all files
        t1 = dspr.variables['time'][:]
        t2 = np.array([dt.timedelta(hours = t)+dt.datetime(1800, 1, 1, 0, 0) for t in t1])    
        #get precipitation rate and close
        #this is in kg/m**2/s (or mm/s)
        time = np.append(time, t2)
        
        #get precipitation rate and close
        #this is in kg/m**2/s (or mm/s)
        prate = np.concatenate((prate, dspr.variables['prate'][:]), axis = 0)
        #this is in kg/m**2
        prwt = np.concatenate((prwt, dspw.variables['pr_wtr'][:]), axis = 0)
        tmp = np.concatenate((tmp, tmpr.variables['air'][:]), axis = 0)
        tmpr.close()            
        dspw.close()
        dspr.close()   
    
    #calculate residence time for each day
    #RT = np.divide(prwt, prate*3600*24) 
    P = np.empty((len(yrs)*12, len(lat), len(lon)))
    monthlytime = []
    Q = copy.copy(P)
    T = copy.copy(P)
    prate = np.ma.array(prate, mask = prate == prate.max())
    prwt = np.ma.array(prwt, mask = prwt == prate.max())
    tmp = np.ma.array(tmp, mask = tmp == prate.max())
    #get the indices for the season
    mons = np.array([m.month for m in time])
    years = np.array([y.year for y in time])
    for i, y in zip(range(len(yrs)), yrs):
        for m in range(1,13):
            inds = np.where((years == y)&(mons == m))[0]
            P[(((i*12)+m)-1), :, :] = np.ma.mean(prate[inds, :, :], axis = 0)
            print(P[(((i*12)+m)-1), :, :])
            Q[(((i*12)+m)-1), :, :] = np.ma.mean(prwt[inds, :, :], axis = 0)
            T[(((i*12)+m)-1), :, :] = np.ma.mean(tmp[inds, :, :], axis = 0)
            monthlytime.append(dt.datetime(y, m, 1, 0, 0))
            
    RT = np.divide(Q, P*3600*24) #makes units days
    #now get seasonal average
    RTseas = np.empty((len(yrs), len(lat), len(lon), 2))
    Tseas = np.empty((len(yrs), len(lat), len(lon), 2))
    mons = np.array([m.month for m in monthlytime])
    years = np.array([y.year for y in monthlytime])
    for i, y in zip(range(len(yrs)), yrs): 
        for sea in seasdict.keys():
            ind = np.where((years == y)&(mons == seasdict[sea][0]))[0]
            subinds = np.array([ind-1, ind, ind+1])
            subinds = subinds[(subinds >= 0)&(subinds < len(RT))]
            RTseas[i, :, :, seasdict[sea][1]] = np.nanmean(RT[subinds, :, :], axis = 0)
            Tseas[i, :, :, seasdict[sea][1]] = np.nanmean(T[subinds, :, :], axis = 0)
            

    return RTseas, Tseas
        

