# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:03:38 2019

@author: u0929173
"""
import numpy as np
import matplotlib as mpl
import os
import cartopy as ctp
import matplotlib.patches as patches
import cartopy.crs as ccrs
from scipy.stats import norm
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

def changemap(glons, glats, regresults, cgseas, cgtcov, cgnmons, years, aggplots):
    RdBu = mpl.pyplot.get_cmap('RdBu_r')
    seasnames = {'JJA':(1, 0), 'DJF':(0, 1)}
    os.chdir(aggplots)  
    
    fig = mpl.pyplot.figure(figsize = (9, 7))
    ax1 = mpl.pyplot.subplot(2, 1, 1, projection=ccrs.Robinson(central_longitude = 0))
    ax2 = mpl.pyplot.subplot(2, 1, 2, projection=ccrs.Robinson(central_longitude = 0))
    axes = np.array([ax1, ax2])
    Lon_e = np.arange(-180, 185, 5)
    Lat_e = np.arange(-90, 95, 5)
    x, y = np.meshgrid(Lon_e, Lat_e)  
    for sn in seasnames.keys():
        ax = axes[seasnames[sn][1]]
        k = seasnames[sn][0]
        #make a plot of all grids and where there is data
        ax.add_feature(ctp.feature.OCEAN, zorder=0, facecolor = '#C2DBE5')
        ax.add_feature(ctp.feature.LAND, zorder=0, edgecolor='dimgray', linewidth = 0.35, facecolor = 'gray')
 
        for i in range(len(glons)):
            ax.plot(x[:, i], y[:, i], color = 'dimgray', linewidth = 0.1, 
                    zorder = 10, transform = ccrs.PlateCarree())    
        for j in range(len(glats)):
            ax.plot(x[j, :], y[j, :], color = 'dimgray', linewidth = 0.1,
                    zorder = 10, transform = ccrs.PlateCarree())
            
        mask = ~np.isnan(regresults[:, :, 4, k])&(regresults[:, :, 4, k] <= 0.1)
        vals = np.ma.array(regresults[:, :, 2, k], mask = mask)
        p = ax.pcolor(x, y, vals, cmap = RdBu, vmin = -1.0, vmax = 1.0, 
                      edgecolors = 'dimgray', linewidths = 0.5, transform = ccrs.PlateCarree())
        mask = ~np.isnan(regresults[:, :, 4, k])&(regresults[:, :, 4, k] > 0.1)
        vals = np.ma.array(regresults[:, :, 2, k], mask = mask)
        p = ax.pcolor(x, y, vals, cmap = RdBu, vmin = -1.0, vmax = 1.0, 
                      edgecolors = 'k', linewidths = 0.8, transform = ccrs.PlateCarree())
        

    axes[0].annotate('(a) JJA', xy = (0.02, 0.94), xytext = (0.05, 0.95), xycoords = 'axes fraction')
    axes[1].annotate('(b) DJF', xy = (0.02, 0.94), xytext = (0.05, 0.95), xycoords = 'axes fraction')
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    cb = fig.colorbar(p, cax=cbar_ax)
    #cb = mpl.pyplot.colorbar(p, orientation = 'horizontal', cax = ax)
    cb.set_label('$\delta^{18}O$ slope ('+u'\u2030'+' decade $^{-1}$)')  
    mpl.pyplot.savefig('Fig1.png', dpi = 300)  
    mpl.pyplot.close()  
    
def PlotModRT(aggplots, Lons, Lats, RTslope, MODslope):
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
    for (ax, k) in zip(axes, [1, 0]):
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
    
    mpl.pyplot.savefig(os.path.join(aggplots, 'Fig4.png'), dpi = 500)
    
def maritime(df, figloc):
    #declare the colormap
    RdGy_r = mpl.pyplot.get_cmap('RdGy_r')
    #get subsets of dataframe
    maritime = df[(df['Maritime'] == 1)&(abs(df['Lat'])<= 15)]
    subset = df[(df['Maritime'] == 1)&(abs(df['Lat'])> 35)]
    #index by season
    winterinds = subset[subset['Season'] == 'DJF'].index
    summerinds = subset[subset['Season'] == 'JJA'].index
    
    fig, axes = mpl.pyplot.subplots(figsize = (7, 4), ncols = 2, sharey = True)
    p = axes[0].scatter(subset.loc[winterinds, 'Press'].values*10, subset.loc[winterinds, 'Slope'].values, 
                       c = subset.loc[winterinds, 'Lat'].values, vmax = 90, vmin = -90, cmap = RdGy_r, 
                       edgecolor = 'k', marker = 's')
    axes[0].errorbar(subset.loc[winterinds, 'Press'].values*10, subset.loc[winterinds, 'Slope'].values, 
                       xerr = subset.loc[winterinds, 'Press_std'].values*10, yerr = subset.loc[winterinds, 'Slope_std'].values, 
                       fmt = '.', ecolor = 'silver', zorder = -1)
    p = axes[0].scatter(subset.loc[summerinds, 'Press'].values*10, subset.loc[summerinds, 'Slope'].values, 
                       c = subset.loc[summerinds, 'Lat'].values, vmax = 90, vmin = -90, cmap = RdGy_r, 
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
                       c = maritime.loc[winterinds, 'Lat'], marker = 's',
                       cmap = RdGy_r, edgecolor = 'k', vmin = -90, vmax = 90)
    axes[1].errorbar(maritime.loc[winterinds, 'Omega']*10, maritime.loc[winterinds, 'Slope'], 
                       xerr = maritime.loc[winterinds, 'Omega_std']*10, yerr = maritime.loc[winterinds, 'Slope_std'], 
                       fmt = '.', ecolor = 'silver', zorder = -1)
    p2 = axes[1].scatter(maritime.loc[summerinds, 'Omega']*10, maritime.loc[summerinds, 'Slope'], 
                       c = maritime.loc[summerinds, 'Lat'], 
                       cmap = RdGy_r, edgecolor = 'k', vmin = -90, vmax = 90)
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
    cb.set_label(r'Latitude $^{\circ}$N')
    
    axes[1].set_xlabel(r'Slope of vertical velocity (Pa s$^{-1}$ decade $^{-1}$)')
    fig.subplots_adjust(top = 0.95, bottom = 0.32, left = 0.12, right = 0.95)
    axes[0].annotate('(a)', xy = (0.02, 0.85), xytext = (0.03, 0.9), xycoords = 'axes fraction')
    axes[1].annotate('(b)', xy = (0.02, 0.85), xytext = (0.03, 0.9), xycoords = 'axes fraction')
    
    mpl.pyplot.savefig(os.path.join(figloc, 'Fig3.png'), dpi = 500)

def stackCompositeFig(hgts, hgtlat, hgtlon, years, stacksNW, stacksNW_std, stacksEur, 
                      stacksEur_std, figloc):
    mask = ~np.isnan(stacksNW)
    corrgridNW = np.zeros((hgts.shape[1], hgts.shape[2]))
    for i in range(len(hgtlat)):
        for j in range(len(hgtlon)):        
            corrgridNW[i, j] = np.corrcoef(hgts[:, i, j][mask], stacksNW[mask])[0,1]
            
    mask = ~np.isnan(stacksEur)
    corrgridEur = np.zeros((hgts.shape[1], hgts.shape[2]))
    for i in range(len(hgtlat)):
        for j in range(len(hgtlon)):        
            corrgridEur[i, j] = np.corrcoef(hgts[:, i, j][mask], stacksEur[mask])[0,1]
            
    levels = [-0.5, 0.5]
    #gs_kw = dict(height_ratios=heights)    
    RdGy = mpl.pyplot.get_cmap('RdGy_r')  
    fig = mpl.pyplot.figure(figsize = (7.5, 4), constrained_layout=True)
    gs = GridSpec(2, 2, height_ratios=[2.5, 1], wspace=0.1, hspace=0.025) 

    #plot the correlation in N hemisphere -- need significance bands too
    hgtlat_edge = np.concatenate(([hgtlat[0]+(2.5/2)], hgtlat-(2.5/2))) 
    x,y = np.meshgrid(np.concatenate((hgtlon, [360])), hgtlat_edge)
    xctr, yctr = np.meshgrid(hgtlon+(2.5/2), hgtlat)

    ax1 = mpl.pyplot.subplot(gs[0], projection=ccrs.PlateCarree())
    ax1.add_feature(ctp.feature.COASTLINE, zorder=10, edgecolor='dimgray', linewidth = 0.75)
    p = ax1.pcolormesh(x, y, corrgridNW,  vmin = -1, vmax = 1, cmap = RdGy, 
                      transform = ccrs.PlateCarree(), alpha = 1)
    cs= ax1.contour(xctr, yctr, corrgridNW, levels, colors = 'gray',
              transform = ccrs.PlateCarree(), linewidth=1)
    ax1.clabel(cs, inline=1, fontsize=10)
    ax1.set_extent([-170, -55, 10, 80], ccrs.PlateCarree())
    
    ax2 = mpl.pyplot.subplot(gs[1], projection=ccrs.PlateCarree())
    ax2.add_feature(ctp.feature.COASTLINE, zorder=10, edgecolor='dimgray', linewidth = 0.75)
    p = ax2.pcolormesh(x, y, corrgridEur,  vmin = -1, vmax = 1, cmap = RdGy, 
                      transform = ccrs.PlateCarree(), alpha = 1)
    cs = ax2.contour(xctr, yctr, corrgridEur, levels, colors = 'gray',
                  transform = ccrs.PlateCarree(), alpha = 1, linewidth=1)
    ax2.clabel(cs, inline=1, fontsize=10)
    ax2.set_extent([-50, 65, 10, 80], ccrs.PlateCarree())

    ax3 = mpl.pyplot.subplot(gs[2])
    
    ax3.errorbar(np.unique(years), stacksNW, yerr = stacksNW_std, fmt = 'o-', 
                        ecolor = '#bababa', markerfacecolor = '#878787', 
                        markeredgecolor = 'black', color = 'black')
    ax3.plot([1960, 2020], [0, 0], color = '#bababa', zorder = 0)
             
    ax4 = mpl.pyplot.subplot(gs[3])         
    ax4.errorbar(np.unique(years), stacksEur, yerr = stacksEur_std, fmt = 'o-', 
                        ecolor = '#bababa', markerfacecolor = '#878787', 
                        markeredgecolor = 'black', color = 'black')
    ax4.plot([1960, 2020], [0, 0], color = '#bababa', zorder = 0)
    ax3.set_ylabel( '$\delta^{18}O$ anomaly'+'('+u'\u2030'+')')
    ax3.annotate('Year', xy = (0.8, 0.02), xytext = (1, -0.4), xycoords = 'axes fraction')
    ax3.set_xlim([1960, 2020])
    ax4.set_xlim([1960, 2020])
    ax3.set_ylim([-2.5, 3])
    ax4.set_ylim([-2.5, 3])
    ax3.spines['top'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax1.annotate('(a)', xy=(0.03, 0.9), xytext = (0.03, 0.9), xycoords = 'axes fraction')
    ax2.annotate('(b)', xy=(0.03, 0.9), xytext = (0.03, 0.9), xycoords = 'axes fraction')
    ax3.annotate('(c)', xy=(0.03, 0.9), xytext = (0.03, 0.9), xycoords = 'axes fraction')
    ax4.annotate('(d)', xy=(0.03, 0.9), xytext = (0.03, 0.9), xycoords = 'axes fraction')
    ax4.set_yticks([])
    ax3.set_xticklabels(ax3.get_xticks().astype(int), rotation = 30)
    ax4.set_xticklabels(ax4.get_xticks().astype(int), rotation = 30)
    
    #make space for a colorbar by  using subplots_adjust
    mpl.pyplot.subplots_adjust(bottom=0.12, right=0.8, top=1.1, left = 0.08)
    cax = mpl.pyplot.axes([0.82, 0.2, 0.025, 0.6])
    cb = mpl.pyplot.colorbar(p, cax=cax)  
    cb.set_label('Correlation coefficient of\n $\delta^{18}O$'+
                 ' anomaly with 500 hPa geopotential height')
    mpl.pyplot.savefig(os.path.join(figloc, 'Fig4.eps'), dpi = 500)
    

  