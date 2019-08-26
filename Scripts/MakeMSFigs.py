# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:03:38 2019

@author: u0929173
"""

def changemap(os, mpl, np, patches, ctp, ccrs, glons, glats, 
                   regresults, cgseas, cgtcov, cgnmons, years, aggplots):
    RdBu = mpl.pyplot.get_cmap('RdBu_r')
    seasnames = ['JJA', 'DJF']
    os.chdir(aggplots)  
    
    fig = mpl.pyplot.figure(figsize = (9, 7))
    ax1 = mpl.pyplot.subplot(2, 1, 1, projection=ccrs.Robinson(central_longitude = 0))
    ax2 = mpl.pyplot.subplot(2, 1, 2, projection=ccrs.Robinson(central_longitude = 0))
    axes = np.array([ax1, ax2])
    Lon_e = np.arange(-180, 185, 5)
    Lat_e = np.arange(-90, 95, 5)
    x, y = np.meshgrid(Lon_e, Lat_e)  
    for (ax, sn, k) in zip(axes, seasnames, range(len(axes))):
        #make a plot of all grids and where there is data
        ax.add_feature(ctp.feature.OCEAN, zorder=0, facecolor = '#C2DBE5')
        ax.add_feature(ctp.feature.LAND, zorder=0, edgecolor='dimgray', linewidth = 0.35, facecolor = 'gray')
 
        for i in range(len(glons)):
            ax.plot(x[:, i], y[:, i], color = 'dimgray', linewidth = 0.1, 
                    zorder = 10, transform = ccrs.PlateCarree())    
        for j in range(len(glats)):
            ax.plot(x[j, :], y[j, :], color = 'dimgray', linewidth = 0.1,
                    zorder = 10, transform = ccrs.PlateCarree())
        vals = np.ma.array(regresults[:, :, 2, k], mask = np.isnan(regresults[:, :, 2, k]))
        p = ax.pcolor(x, y, vals, cmap = RdBu, vmin = -1.0, vmax = 1.0, 
                      edgecolors = 'k', linewidths = 0.5, transform = ccrs.PlateCarree())

    axes[0].annotate('(a) '+seasnames[0], xy = (0.02, 0.94), xytext = (0.05, 0.95), xycoords = 'axes fraction')
    axes[1].annotate('(b) '+seasnames[1], xy = (0.02, 0.94), xytext = (0.05, 0.95), xycoords = 'axes fraction')
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    cb = fig.colorbar(p, cax=cbar_ax)
    #cb = mpl.pyplot.colorbar(p, orientation = 'horizontal', cax = ax)
    cb.set_label('$\delta^{18}O$ slope ('+u'\u2030'+' decade $^{-1}$)')  
    mpl.pyplot.savefig('Fig1.png', dpi = 300)  
    mpl.pyplot.close()  
    
def PlotModRT(os, np, mpl, ctp, ccrs, aggplots, Lons, Lats, RTslope, MODslope):
    RdBu = mpl.pyplot.get_cmap('RdBu_r')
    Lon_e = np.arange(-180, 185, 5)
    Lat_e = np.arange(-90, 95, 5)
    x, y = np.meshgrid(Lon_e, Lat_e)    
    ss = (((RTslope > 0)*(MODslope > 0))+((RTslope < 0)*(MODslope < 0))).astype(int)
    ss = np.ma.masked_less(ss, 0.5)

    #stipple where the sign is the same
    fig = mpl.pyplot.figure(figsize = (9, 7))
    ax1 = mpl.pyplot.subplot(2, 1, 1, projection=ccrs.Robinson(central_longitude = 0))
    ax2 = mpl.pyplot.subplot(2, 1, 2, projection=ccrs.Robinson(central_longitude = 0))
    axes = np.array([ax1, ax2])
    for (ax, k) in zip(axes, range(len(axes))):
        #make a plot of all grids and where there is data
        ax.add_feature(ctp.feature.COASTLINE, zorder=10)
 
        vals = np.ma.array(MODslope[:, :, k], mask = np.isnan(MODslope[:, :, k]))
        p = ax.pcolormesh(x, y, vals, cmap = RdBu, vmin = -1, vmax = 1, 
                      transform = ccrs.PlateCarree())  
        ax.pcolor(x, y, ss[:, :, k], hatch = '////', zorder = 5, alpha = 0.0, transform = ccrs.PlateCarree())

    axes[0].annotate('(a) JJA', xy = (0.02, 0.95), xytext = (0.02, 0.95), xycoords = 'axes fraction')    
    axes[1].annotate('(b) DJF', xy = (0.02, 0.95), xytext = (0.02, 0.95), xycoords = 'axes fraction')
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    cb = fig.colorbar(p, cax=cbar_ax)
    #cb = mpl.pyplot.colorbar(p, orientation = 'horizontal', cax = ax)
    cb.set_label('$\delta^{18}O$ slope ('+u'\u2030'+' decade $^{-1}$)')
    
    mpl.pyplot.savefig(os.path.join(aggplots, 'FigS1.png'), dpi = 500)
    
def StackedBar(os, np, mpl, patches, df, figloc):
    df['sig'] = np.multiply(np.add(df['Slope'], df['Slope_std']),np.add(df['Slope'], -df['Slope_std'])) > 0

    df['Colors'] = np.where(df['Maritime'].values == 1, '#35978f', '#bf812d')
      #change color of significant sites
    df['Colors'] = np.where((df['Colors'] == '#35978f')&(df['sig']), '#01665e', df['Colors'])
    df['Colors'] = np.where((df['Colors'] == '#bf812d')&(df['sig']), '#8c510a', df['Colors'])
    
    latbins = np.array([-66, -43, -23, 0, 23, 43, 66])
    #for plotting
    latctrs_c = np.array([-52, -29, -8, 16, 36, 59])
    latctrs_m = np.array([-59, -36, -16, 8, 29, 52])
    fig, axes = mpl.pyplot.subplots(nrows = 2, ncols = 1, sharex = True, figsize = (7, 5.5))
    default = 0.01
    width = 5
    edgecolors = ['dimgray', 'black']
    for axcts in range(2):
        if axcts == 0:
            key = 'JJA'
        else:
            key = 'DJF'
        for i in range(len(latctrs_c)):
            latbininds = df.index[(df['Lat'].values >= latbins[i]) & (df['Lat'].values <= latbins[i+1])]
            for j in range(2):
                subinds = latbininds[(df['Maritime'][latbininds] == 0)&(df['Season'][latbininds] == key)&(df['sig'][latbininds] == bool(j))]
                subinds_n = subinds[df['Slope'][subinds] < 0]
                subinds_p = subinds[df['Slope'][subinds] > 0]
                if len(subinds_n) > 0:
                    if len(subinds_n) > 1:
                        position = latctrs_c[i]
                        llc = (min(df['Slope'][subinds_n].values), position)
                        hgt = max(df['Slope'][subinds_n].values) - min(df['Slope'][subinds_n].values)
                        if hgt < default:
                                hgt = default
                        color = np.unique(df['Colors'][subinds_n].values)[0]
                        pat = patches.Rectangle(xy = llc, width = hgt, height = width, 
                                fill = True, facecolor = color, edgecolor = edgecolors[j], zorder = j)
                        axes[axcts].add_patch(pat) 
                    else:
                        axes[axcts].errorbar(df.loc[subinds_n, 'Slope'], latctrs_c[i]+(width/2), 
                            xerr = df.loc[subinds_n, 'Slope_std'], fmt = 'o', ecolor = edgecolors[j], 
                            zorder = 0)
                        axes[axcts].scatter(df.loc[subinds_n, 'Slope'], latctrs_c[i]+(width/2), 
                            c = df.loc[subinds_n, 'Colors'], s = 40, edgecolor = edgecolors[j], zorder = 3)
                if len(subinds_p) > 0:
                    if len(subinds_p) > 1:
                        position = latctrs_c[i]
                        llc = (min(df['Slope'][subinds_p].values), position)
                        hgt = max(df['Slope'][subinds_p].values) - min(df['Slope'][subinds_p].values)
                        if hgt < default:
                                hgt = default
                        color = np.unique(df['Colors'][subinds_p].values)[0]
                        pat = patches.Rectangle(xy = llc, width = hgt, height = width, 
                                fill = True, facecolor = color, edgecolor = edgecolors[j], zorder = j)
                        axes[axcts].add_patch(pat) 
                    else:
                        axes[axcts].errorbar(df.loc[subinds_p, 'Slope'], latctrs_c[i]+(width/2), 
                            xerr = df.loc[subinds_p, 'Slope_std'], fmt = 'o', ecolor = edgecolors[j], 
                            zorder = 0)
                        axes[axcts].scatter(df.loc[subinds_p, 'Slope'], latctrs_c[i]+(width/2), 
                            c = df.loc[subinds_p, 'Colors'], s = 40, edgecolor = edgecolors[j], zorder = 3)
    
    
                subinds = latbininds[(df['Maritime'][latbininds] == 1)&(df['Season'][latbininds] == key)&(df['sig'][latbininds] == bool(j))]
                subinds_n = subinds[df['Slope'][subinds] < 0]
                subinds_p = subinds[df['Slope'][subinds] > 0]
                if len(subinds_n) > 0:
                    if len(subinds_n) > 1:
                        position = latctrs_m[i]
                        llc = (min(df['Slope'][subinds_n].values), position)
                        hgt = max(df['Slope'][subinds_n].values) - min(df['Slope'][subinds_n].values)
                        if hgt < default:
                                hgt = default
                        color = np.unique(df['Colors'][subinds_n].values)[0]
                        pat = patches.Rectangle(xy = llc, width = hgt, height = width, 
                                fill = True, facecolor = color, edgecolor = edgecolors[j], zorder = j)
                        axes[axcts].add_patch(pat) 
                    else:
                        axes[axcts].errorbar(df.loc[subinds_n, 'Slope'], latctrs_m[i]+(width/2), 
                            xerr = df.loc[subinds_n, 'Slope_std'], fmt = 'o', ecolor = edgecolors[j], 
                            zorder = 0)
                        axes[axcts].scatter(df.loc[subinds_n, 'Slope'], latctrs_m[i]+(width/2), 
                            c = df.loc[subinds_n, 'Colors'], s = 40, edgecolor = edgecolors[j], zorder = 3)
                if len(subinds_p) > 0:
                    if len(subinds_p) > 1:
                        position = latctrs_m[i]
                        llc = (min(df['Slope'][subinds_p].values), position)
                        hgt = max(df['Slope'][subinds_p].values) - min(df['Slope'][subinds_p].values)
                        if hgt < default:
                                hgt = default
                        color = np.unique(df['Colors'][subinds_p].values)[0]
                        pat = patches.Rectangle(xy = llc, width = hgt, height = width, 
                                fill = True, facecolor = color, edgecolor = edgecolors[j], zorder = j)
                        axes[axcts].add_patch(pat)  
                    else:
                        axes[axcts].errorbar(df.loc[subinds_p, 'Slope'], latctrs_m[i]+(width/2), 
                            xerr = df.loc[subinds_p, 'Slope_std'], fmt = '.', ecolor = edgecolors[j], 
                            zorder = 0)
                        axes[axcts].scatter(df.loc[subinds_p, 'Slope'], latctrs_m[i]+(width/2), 
                            c = df.loc[subinds_p, 'Colors'], s = 40, edgecolor = edgecolors[j], zorder = 3)
                    
    inds = df.index[(df['Lat'].values >= 66)&(df['Season'] == 'JJA')]
    axes[0].errorbar(df.loc[inds, 'Slope'], df.loc[inds, 'Lat'], 
                            xerr = df.loc[inds, 'Slope_std'], fmt = '.', ecolor = 'black', 
                            zorder = 0)
    axes[0].scatter(df.loc[inds, 'Slope'], df.loc[inds, 'Lat'], c = df.loc[inds, 'Colors'], s = 40, edgecolor = 'k')
    inds = df.index[(df['Lat'].values <= -66)&(df['Season'] == 'JJA')]
    axes[0].errorbar(df.loc[inds, 'Slope'], df.loc[inds, 'Lat'], 
                            xerr = df.loc[inds, 'Slope_std'], fmt = '.', ecolor = 'black', 
                            zorder = 0)
    axes[0].scatter(df.loc[inds, 'Slope'], df.loc[inds, 'Lat'], c = df.loc[inds, 'Colors'], s = 40, edgecolor = 'k')
    
    inds = df.index[(df['Lat'].values >= 66)&(df['Season'] == 'DJF')]
    axes[1].errorbar(df.loc[inds, 'Slope'], df.loc[inds, 'Lat'], 
                            xerr = df.loc[inds, 'Slope_std'], fmt = '.', ecolor = 'black', 
                            zorder = 0)
    axes[1].scatter(df.loc[inds, 'Slope'], df.loc[inds, 'Lat'], c = df.loc[inds, 'Colors'], s = 40, edgecolor = 'k')
    inds = df.index[(df['Lat'].values <= -66)&(df['Season'] == 'DJF')]
    axes[1].errorbar(df.loc[inds, 'Slope'], df.loc[inds, 'Lat'], 
                            xerr = df.loc[inds, 'Slope_std'], fmt = '.', ecolor = 'black', 
                            zorder = 0)
    axes[1].scatter(df.loc[inds, 'Slope'], df.loc[inds, 'Lat'], c = df.loc[inds, 'Colors'], s = 40, edgecolor = 'k')
    
    axes[0].plot([0,0], [-90, 90], linestyle = '--', color = 'silver', linewidth = 1.5, zorder = -2)
    axes[1].plot([0,0], [-90, 90], linestyle = '--', color = 'silver', linewidth = 1.5, zorder = -2)
    axes[0].set_xlim([-2, 2])
    axes[1].set_ylim([-90, 90])
    axes[0].set_ylim([-90, 90])
    axes[1].set_xlabel('$\delta^{18}O$ slope ('+u'\u2030'+' decade $^{-1}$)')
    axes[1].annotate('Latitude', xy = (0.02, 0.5), xytext = (0.04, 0.52), xycoords = 'figure fraction', rotation = 'vertical')
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[0].spines['top'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[0].annotate('(a) JJA', xy = (0.02, 0.95), xytext = (0.03, 0.9), xycoords = 'axes fraction')
    axes[1].annotate('(b) DJF', xy = (0.02, 0.95), xytext = (0.03, 0.9), xycoords = 'axes fraction')
    
    markers = [mpl.pyplot.Line2D([0,0],[0,0], color = col, marker='o', linestyle='') 
        for col in ['#01665e', '#35978f', '#8c510a', '#bf812d']]
    axes[0].legend(np.array(markers), np.array(['Maritime: significant', 'Maritime: not significant', 
        'Continental: significant', 'Continental: not significant']), numpoints=1, loc = 'upper left', 
            ncol = 2, fontsize = 10, bbox_to_anchor=(0.05, 1.3), frameon = False)
    mpl.pyplot.savefig(os.path.join(figloc, 'Fig2.png'), dpi = 500)

def maritime(os, sm, mpl, df, figloc):
    #declare the colormap
    RdGy_r = mpl.pyplot.get_cmap('RdGy_r')
    #get subsets of dataframe
    maritime = df[(df['Maritime'] == 1)&(abs(df['Lat'])<= 20)]
    subset = df[(df['Maritime'] == 1)&(abs(df['Lat'])> 25)]
    #index by season
    winterinds = subset[subset['Season'] == 'DJF'].index
    summerinds = subset[subset['Season'] == 'JJA'].index
    
    fig, axes = mpl.pyplot.subplots(figsize = (7, 4), ncols = 2, sharey = True)
    p = axes[0].scatter(subset.loc[winterinds, 'Press'].values*10, subset.loc[winterinds, 'Slope'].values, 
                       c = subset.loc[winterinds, 'T'].values*10, vmax = 0.5, vmin = -0.5, cmap = RdGy_r, 
                       edgecolor = 'k', marker = 's')
    axes[0].errorbar(subset.loc[winterinds, 'Press'].values*10, subset.loc[winterinds, 'Slope'].values, 
                       xerr = subset.loc[winterinds, 'Press_std'].values*10, yerr = subset.loc[winterinds, 'Slope_std'].values, 
                       fmt = '.', ecolor = 'silver', zorder = -1)
    p = axes[0].scatter(subset.loc[summerinds, 'Press'].values*10, subset.loc[summerinds, 'Slope'].values, 
                       c = subset.loc[summerinds, 'T'].values*10, vmax = 0.5, vmin = -0.5, cmap = RdGy_r, 
                       edgecolor = 'k')
    axes[0].errorbar(subset.loc[summerinds, 'Press'].values*10, subset.loc[summerinds, 'Slope'].values, 
                       xerr = subset.loc[summerinds, 'Press_std'].values*10, yerr = subset.loc[summerinds, 'Slope_std'].values, 
                       fmt = '.', ecolor = 'silver', zorder = -1)

    axes[0].set_xlabel(r'Slope of Surface Pressure (Pa decade $^{-1}$)')
    axes[0].set_ylabel(r'Slope of $\delta^{18}O$ ('+u'\u2030'+' decade $^{-1}$)')
    axes[0].set_ylim([-1.5, 1.15])
    axes[0].set_xlim([-110, 100])
    axes[0].plot([-125, 100], [0, 0], color = 'dimgray', zorder = -2, linewidth = 2, linestyle = '--')
    axes[0].plot([0, 0], [-1.5, 1.15], color = 'dimgray', zorder = -2, linewidth = 2, linestyle = '--')
    
    winterinds = maritime[maritime['Season'] == 'DJF'].index
    summerinds = maritime[maritime['Season'] == 'JJA'].index
    
    p1 = axes[1].scatter(maritime.loc[winterinds, 'Omega']*10, maritime.loc[winterinds, 'Slope'], 
                       c = maritime.loc[winterinds, 'T']*10, marker = 's',
                       cmap = RdGy_r, edgecolor = 'k', vmin = -0.5, vmax = 0.5)
    axes[1].errorbar(maritime.loc[winterinds, 'Omega']*10, maritime.loc[winterinds, 'Slope'], 
                       xerr = maritime.loc[winterinds, 'Omega_std']*10, yerr = maritime.loc[winterinds, 'Slope_std'], 
                       fmt = '.', ecolor = 'silver', zorder = -1)
    p2 = axes[1].scatter(maritime.loc[summerinds, 'Omega']*10, maritime.loc[summerinds, 'Slope'], 
                       c = maritime.loc[summerinds, 'T']*10, 
                       cmap = RdGy_r, edgecolor = 'k', vmin = -0.5, vmax = 0.5)
    axes[1].errorbar(maritime.loc[summerinds, 'Omega']*10, maritime.loc[summerinds, 'Slope'], 
                       xerr = maritime.loc[summerinds, 'Omega_std']*10, yerr = maritime.loc[summerinds, 'Slope_std'], 
                       fmt = '.', ecolor = 'silver', zorder = -1)
    axes[1].plot([-0.01, 0.01], [0, 0], color = 'dimgray', zorder = -2, linewidth = 2, linestyle = '--')
    axes[1].plot([0, 0], [-1.5, 1.15], color = 'dimgray', zorder = -2, linewidth = 2, linestyle = '--')
    axes[1].legend([p1, p2], ['DJF', 'JJA'], ncol = 2, loc = (0.38, 0.05))
    
    lg = axes[1].get_legend()
    lg.legendHandles[1].set_facecolor('k')
    lg.legendHandles[0].set_facecolor('k')
    axes[1].set_xlim([-0.01, 0.01])
    axes[1].set_ylim([-1.5, 1.15])
    
    cbar_ax2 = fig.add_axes([0.25, 0.12, 0.50, 0.03])
    cb = fig.colorbar(p1, cax=cbar_ax2, orientation = 'horizontal')
    cb.set_label(r'Slope  of T ($^{\circ}$ C decade $^{-1}$)')
    
    axes[1].set_xlabel(r'Slope of vertical velocity (Pa s$^{-1}$ decade $^{-1}$)')
    fig.subplots_adjust(top = 0.95, bottom = 0.32, left = 0.12, right = 0.95)
    axes[0].annotate('(a)', xy = (0.02, 0.85), xytext = (0.03, 0.9), xycoords = 'axes fraction')
    axes[1].annotate('(b)', xy = (0.02, 0.85), xytext = (0.03, 0.9), xycoords = 'axes fraction')
    
    mpl.pyplot.savefig(os.path.join(figloc, 'Fig3.png'), dpi = 500)

def MakeGeopotentialRegress(os, np, mpl, ctp, ccrs, stack, hgts, hgtlat, hgtlon, stackname, hgtname, figloc):
    #correlation plot of 500 mbar geopotential height with NW stack- this one is pretty convincing
    mask = ~np.isnan(stack)
    corrgrid = np.zeros((hgts.shape[1], hgts.shape[2]))
    for i in range(len(hgtlat)):
        for j in range(len(hgtlon)):        
            corrgrid[i, j] = np.corrcoef(hgts[:, i, j][mask], stack[mask])[0,1]
    
    RdGy = mpl.pyplot.get_cmap('RdGy_r')
    #plot the correlation in N hemisphere -- need significance bands too
    x,y = np.meshgrid(hgtlon, hgtlat)
    fig = mpl.pyplot.figure(figsize = (7, 5.5))
    ax = mpl.pyplot.axes(projection = ccrs.Robinson(central_longitude = -100))
    ax.add_feature(ctp.feature.COASTLINE, zorder=10, edgecolor='dimgray', linewidth = 0.5)
    p = ax.pcolormesh(x, y, corrgrid,  vmin = -1, vmax = 1, cmap = RdGy, transform = ccrs.PlateCarree(), alpha = 1)
    cb = mpl.pyplot.colorbar(p)  
    cb.set_label('Correlation coefficient') 
    ax.set_title('{} stack vs \n{} mbar geopotential height'.format(stackname, hgtname))   
    mpl.pyplot.savefig(os.path.join(figloc, '{}_{}hgt.png'.format(stackname, hgtname)), dpi = 500)  
    mpl.pyplot.close() 

def MakeStackFig(os, np, mpl, years, stackNW, stackNW_std, stackSE, stackSE_std, PNA_DJF, figloc):
    fig, ax = mpl.pyplot.subplots(nrows = 3, ncols = 1, sharex = True)
    ax[0].errorbar(np.unique(years), stackNW, yerr = stackNW_std, fmt = 'o-', 
                        ecolor = '#bababa', markerfacecolor = '#878787', 
                        markeredgecolor = 'black', color = 'black')
    ax[0].plot([1960, 2016], [0, 0], color = '#bababa', zorder = 0)
    ax[1].errorbar(np.unique(years), stackSE, yerr = stackSE_std, fmt = 'o-', 
                        ecolor = '#f4a582', markerfacecolor = '#d6604d', 
                        markeredgecolor = 'black', color = '#d6604d')
    ax[1].plot([1960, 2016], [0, 0], color = '#bababa', zorder = 0)
    ax[2].plot(np.unique(years), PNA_DJF, marker = 'o', 
                         markerfacecolor = '#67a9cf', 
                        markeredgecolor = 'black', color = '#2166ac')   
    ax[2].plot([1960, 2016], [0, 0], color = '#bababa', zorder = 0)
      
      
    ax[0].annotate('(a) northwest', xy = (0, 100), xytext = (0.02, 0.9), xycoords = 'axes fraction')
    ax[1].annotate('(b) southeast', xy = (0, 100), xytext = (0.02, 0.9), xycoords = 'axes fraction')
    ax[2].annotate('(c) PNA index', xy = (0, 100), xytext = (0.02, 0.9), xycoords = 'axes fraction')
    ax[1].annotate( 'd18O anomaly stack', xy = (-0.05, 0.75), xytext = (0.05, 0.75), 
      xycoords = 'figure fraction', rotation = 'vertical')
    ax[1].set_xlabel('Year')
    ax[2].set_ylabel('PNA index')
    ax[2].set_xlim([1960, 2016])
    mpl.pyplot.savefig(os.path.join(figloc, 'StackPlotPNA.png'), dpi = 500)

def MakeSemivariogramFig(os, np, mpl, semidata, figloc):
    rearth = 6371.0
    timename = {'JJA':0, 'DJF':1} 
    #length scale of synoptic systems
    synoptic = 1000/rearth
    longwave1 = 6000/rearth
    longwave2 = 8000/rearth
    quadrant = np.pi/2

    #'JJA 5x5 agg data', 'JJA 5x5 ens', 'JJA 5x5 ens sub'
    mpl.pyplot.figure(figsize = (4.5, 4))
    xtickdict = {'{} 5x5 agg data': 1, '{} 5x5 ens sub': 2, '{} 5x5 ens': 3}
    for i, sym in zip(timename.keys(), ['o', 's']) :
        for j in xtickdict.keys():
            mpl.pyplot.scatter(xtickdict[j], semidata.loc[1, j.format(i)],
                                    c = '#5ab4ac', edgecolor = 'k', marker = sym, s = 60)
            mpl.pyplot.scatter(xtickdict[j], semidata.loc[0, j.format(i)], 
                                    c = '#d8b365', edgecolor = 'k', marker = sym, s = 60) 
                                    
    p1 = mpl.pyplot.scatter(1, -2, marker = 's', c = '#5ab4ac', edgecolor = 'k', label = 'color: Timeslice, symbol: DJF')  
    p2 = mpl.pyplot.scatter(1, -2, marker = 'o', c = '#d8b365', edgecolor = 'k', label = 'color: Difference, symbol: JJA')
    
    #mpl.pyplot.plot([5.5, 5.5], [0, 2.5], color = 'gray', linewidth = 1.5)   
    mpl.pyplot.plot([0, 7], [synoptic, synoptic], color = 'gray', linewidth = 1.5, zorder = 0, linestyle = '--')
    mpl.pyplot.annotate('Mesoscale', xy = (5, synoptic), xytext = (4.4, synoptic-0.1), xycoords = 'data', size = 8)
    mpl.pyplot.annotate('Synoptic scale', xy = (5, synoptic), xytext = (4.4, synoptic+0.05), xycoords = 'data', size = 8)
    mpl.pyplot.plot([0, 7], [longwave1, longwave1], color = 'gray', linewidth = 1.5, zorder = 0, linestyle = '--')
    mpl.pyplot.annotate('Short waves', xy = (5, longwave1-0.3), xytext = (4.4, longwave1-0.3), xycoords = 'data', size = 8)
    mpl.pyplot.annotate('Rossby waves', xy = (5, longwave2-0.2), xytext = (4.4, longwave2-0.2), xycoords = 'data', size = 8)
    mpl.pyplot.plot([0, 7], [longwave2, longwave2], color = 'gray', linewidth = 1.5, zorder = 0, linestyle = '--')
    mpl.pyplot.plot([0, 7], [quadrant, quadrant], color = 'gray', linewidth = 1.5, zorder = 0, linestyle = '--')
    mpl.pyplot.annotate('Globe quadrant', xy = (5, quadrant-0.05), xytext = (4.4, quadrant-0.1), xycoords = 'data', size = 8)
    #mpl.pyplot.plot([0, 7], [hemisphere, hemisphere], color = 'gray', linewidth = 1.5, zorder = 0)
                         
    mpl.pyplot.ylim([0, 1.75])
    mpl.pyplot.xlim([0, 6])
    #mpl.pyplot.xlabel('dataset')
    labels = [lab[3:] for lab in xtickdict.keys()]
    mpl.pyplot.xticks(xtickdict.values(), labels, rotation = 20)
    mpl.pyplot.ylabel('Semivariogram lag, radians')
    mpl.pyplot.legend([p1, p2], ['color: average, symbol: DJF', 'color: change, symbol: JJA'], 
                      loc = (0.15, 0.93), framealpha = 1)
    mpl.pyplot.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.95)
    mpl.pyplot.savefig(os.path.join(figloc, 'SemivariogramLag.png'), dpi = 500)