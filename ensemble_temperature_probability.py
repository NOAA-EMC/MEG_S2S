# ==========================================
# Title: 2-m Temperature Outlook Plots
# Description: Plot temperature probability outlooks using ensembles grib files and analysis grid values... 
# outlooks defined as probabilities of above or below normal
# Author: Marcel Caron
# Created: August 11 2021
# Last Modified: Feb. 4 2022
# ==========================================

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patheffects as pe
import cartopy.crs as ccrs
import cartopy.feature as cf
from datetime import datetime, timedelta as td
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata
import scipy.stats as stats
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
import seaborn as sns
import pickle

# ============= CHANGES BELOW ==============

DATA_DIR = '/scratch2/NCEPDEV/stmp1/Shannon.Shields/scripts/s2s_data/'
SAVE_DIR = '/scratch2/NCEPDEV/stmp1/Shannon.Shields/scripts/s2s/cases/webinar/'

# ~~~~~~~
# This is the continuation of DATA_DIR, because each model directory is unique
# Use standard C datetime format codes (e.g., %Y for year) for datetimes, 
# and use PP in place of ensemble mem, FF in place of lead time, STAT in place of climo statistic (mean, stdv ...)
# ~~~~~~~
gefs_string = 'gefs.%Y%m%d/gePP.t%Hz.pgrb2a.0p50.fFF'
rap_string = 'rap.%Y%m%d/rap.t%Hz.awip32fFF.grib2' 
climo_string = 'climo_files/STAT_%m%d' 

# ~~~~~~~
# Valid range, inclusive, using CPC definitions for "Day X"
# ~~~~~~~
day_range = (6, 10) #8, 14 

# ~~~~~~~
# e.g., ['gefs'] 
# ~~~~~~~
models = ['gefs'] 

# ~~~~~~~
# Analysis field used for comparison
# ~~~~~~~
analysis = 'rap'

# ~~~~~~~
# Create a plot for each model at the following init times
# ~~~~~~~
inits = [datetime(2021,12,1,00,0,0)]

# ~~~~~~~
# List of domains to plot
# e.g., northamer, conus, econus, zoom, neast, splains_zoom, other
# ~~~~~~~
domains = ['conus'] 

# ~~~~~~~
# if domain = 'other', will use the parameters below to plot, otherwise ignores
# ~~~~~~~
latrange = (26., 44.)
lonrange = (-109., -84.5)
parallels = (25., 35.)
merid = -98.
figsize = (9.5,7.) # tuple width, height for desired dimensions of plots

# ~~~~~~~
# Put grid resolution in degrees, or put None
# ~~~~~~~
finer_grid = 1 #None

# ============= CHANGES ABOVE ==============


# Works generally if you're plotting in the western hemisphere
def get_lons(da):
   if np.max(da.longitude) > 180.:
      lons = np.subtract(np.mod(np.add(np.array(da.longitude),180.),360.),180.)
   else:
      lons = np.array(da.longitude)
   return lons

def get_domains(domain='northamer', parallels=None, merid=None, latrange=None, lonrange=None, num=0, figsize=None, return_box=False):
   if str(domain).upper() == 'NORTHAMER':
      if not return_box:
         fig = plt.figure(num=num, figsize=(7.6,7.))
      parallels = (23., 37.)
      merid = -95.
      latrange = (10., 74.)
      lonrange = (-143., -49.)
   elif str(domain).upper() == 'WASH':
      if not return_box:
         fig = plt.figure(num=num, figsize=(8.3,7.))
      parallels = (25., 35.)
      merid = -120.
      latrange = (44., 49.5)
      lonrange = (-125., -116.)
   elif str(domain).upper() == 'CONUS':
      if not return_box:
         fig = plt.figure(num=num, figsize=(8.9,7.))
      parallels = (25., 35.)
      merid = -97.5
      latrange = (19., 57.)
      lonrange = (-123., -70.)
   elif str(domain).upper() == 'NEAST':
      if not return_box:
         fig = plt.figure(num=num, figsize=(7.6,7.))
      parallels = (25., 35.)
      merid = -77.
      latrange = (35.1, 47.1)
      lonrange = (-85.5, -69.)
   elif str(domain).upper() == 'SPLAINS_ZOOM':
      if not return_box:
         fig = plt.figure(num=num, figsize=(9.5,7.))
      parallels = (25., 35.)
      merid = -98.
      latrange = (30.5, 40.)
      lonrange = (-104.5, -89.)
   elif str(domain).upper() == 'ZOOM':
      if not return_box:
         fig = plt.figure(num=num, figsize=(7.6,7.))
      parallels = (25., 35.)
      merid = -90.
      latrange = (25., 50.)
      lonrange = (-108., -72.)
   elif str(domain).upper() == 'ECONUS':
      if not return_box:
         fig = plt.figure(num=num, figsize=(9.5,7.))
      parallels = (25., 35.)
      merid = -90.
      latrange = (23.6, 47.4)
      lonrange = (-105., -62.)
   elif str(domain).upper() == 'OTHER':
      if not return_box:
         fig = plt.figure(num=num, figsize=figsize)
      parallels = parallels
      merid = merid
      latrange = latrange
      lonrange = lonrange
   proj = ccrs.LambertConformal(central_longitude=merid, standard_parallels=parallels)
   if return_box:
      return lonrange, latrange
   else:
      return proj, fig

def climo_lower_data(init, climo_string, fstart, fend):
   
   print('========================================')
   #fvalids = np.arange(fstart, fend+1, 6)
   fvalids = np.arange(fstart, fend+1, 24)
   valid_hours = [(init + td(hours=int(hrs))).hour for hrs in fvalids]
   data_upper_temps = []
   data_lower_temps = []
   hours_idx = {'18':0, '0':1, '6':2, '12':3}
   for f, fvalid in enumerate(fvalids):
      if valid_hours[f] in [0,6,12]:
         use_date = init+td(hours=int(fvalid))
      else:
         use_date = init+td(hours=int(fvalid)+24) # 18Z data are stored in the file named after the following day for some reason
      fpath_mean = os.path.join(DATA_DIR, use_date.strftime(climo_string.replace('STAT','mean')))
      fpath_stdev = os.path.join(DATA_DIR, use_date.strftime(climo_string.replace('STAT','stdev')))
      with xr.open_dataset(fpath_mean, engine='cfgrib', backend_kwargs=dict(
         indexpath='',
         filter_by_keys={'typeOfLevel':'heightAboveGround','shortName':'2t'})) as ds:
         temp = ds.t2m
         if f == 0:
            lons = get_lons(temp)
            lats = np.array(temp.latitude)
         data_mean_temp = np.array(temp)
      with xr.open_dataset(fpath_stdev, engine='cfgrib', backend_kwargs=dict(
         indexpath='',
         filter_by_keys={'typeOfLevel':'heightAboveGround','shortName':'2t'})) as ds:
         #temp = ds.mx2t6[hours_idx[str(valid_hours[f])]]
         temp = ds.t2m
         data_stdev_temp = np.array(temp)
      multiplier=.43072 

      # 33.333% of values in a normal distribution fall within +/- .43072 standard deviations from the mean 
      # (we're assuming normal) .43 is used in E. Becker and van den Dool, H. (2016) Probabilistic seasonal 
      # forecasts in the North American Multimodel Ensemble: A baseline skill assessment. J. Climate, 29(8), 3015-26
      data_upper_temp = np.add(data_mean_temp, data_stdev_temp*multiplier)
      data_lower_temp = np.subtract(data_mean_temp, data_stdev_temp*multiplier)
      data_upper_temps.append(data_upper_temp)
      data_lower_temps.append(data_lower_temp)
         
   data_upper = np.average(data_upper_temps, axis=0)
   data_lower = np.average(data_lower_temps, axis=0)

   # convert K to F
   data_upper = np.add(np.multiply(np.subtract(data_upper, 273.15), 9./5.), 32.)
   data_lower = np.add(np.multiply(np.subtract(data_lower, 273.15), 9./5.), 32.)
   
   print('Data Retrieved')
   print('========================================') 
   return (data_upper, data_lower, lats, lons)

def ens_lower_data(init, model, model_string, fstart, fend, ens_mem=None, analysis=False):
   
   print('========================================')
   if ens_mem:
      print(f'Getting and Summing {model} {ens_mem} Member Data')
   else:
      print(f'Getting and Summing {model} Data')
   if analysis:
      #fvalids = np.arange(fstart, fend+1, 12)
      fvalids = np.arange(fstart, fend+1, 24)
   else:
      #fvalids = np.arange(fstart, fend+1, 6)
      fvalids = np.arange(fstart, fend+1, 24)
   valid_hours = [(init + td(hours=int(hrs))).hour for hrs in fvalids]

   data_temps = []
   if str(model).upper() in ['GEFS']: # global models have 3 leading zeros
      for f, fvalid in enumerate(fvalids):
         if analysis:
            fpath = os.path.join(DATA_DIR, (init+td(hours=int(fvalid))).strftime(model_string).replace('FF','{0:0=3d}'.format(0)).replace('PP',ens_mem)) 
         else:
            fpath = os.path.join(DATA_DIR, init.strftime(model_string).replace('FF','{0:0=3d}'.format(int(fvalid))).replace('PP',ens_mem)) 
         with xr.open_dataset(fpath, engine='cfgrib', backend_kwargs=dict(
            indexpath='',
            filter_by_keys={'typeOfLevel':'heightAboveGround'})) as ds:
            temp = ds.t2m
            if f == 0:
               lons = get_lons(temp)
               lats = np.array(temp.latitude)
            data_temp = np.array(temp)
            data_temps.append(data_temp)
   else:
      for f, fvalid in enumerate(fvalids):
         #fpath = os.path.join(DATA_DIR, init.strftime(model_string).replace('FF','{0:0=2d}'.format(int(fvalid))))
         fpath = os.path.join(DATA_DIR, (init+td(hours=int(fvalid))).strftime(model_string).replace('FF','{0:0=2d}'.format(0)))
         with xr.open_dataset(fpath, engine='cfgrib', backend_kwargs=dict(
            indexpath='',
            filter_by_keys={'typeOfLevel':'heightAboveGround'})) as ds:
            temp = ds.t2m
            if f == 0:
               lons = get_lons(temp)
               lats = np.array(temp.latitude)
            data_temp = np.array(temp)
            data_temps.append(data_temp)

   data_temps = np.array(data_temps)
   
   data = np.average(data_temps, axis=0) # mean daily temperature
   data = np.add(np.multiply(np.subtract(data, 273.15), 9./5.), 32.) # convert K to F
   
   print('Data Retrieved')
   print('========================================') 
   return (data, lats, lons)

# Returns numpy data meshgrid as a bilinear interpolation of model data onto analysis grid
def diff_data(anl_lons, anl_lats, model_lons, model_lats, model_data):
   anl_lons, anl_lats, model_lons, model_lats, model_data = [np.array(item) for item in [anl_lons, anl_lats, model_lons, model_lats, model_data]]
   model_data_interpd = griddata(
      (model_lons.flatten(), model_lats.flatten()),
      model_data.flatten(),
      (anl_lons, anl_lats), method='linear')
   return model_data_interpd

def BSS(prob, anal):
   climo = np.divide(np.ones(anal.shape), 3) # .33 probability that fcst is in any given tercile (above, below, or normal), for all locations
   BS = np.nanmean(np.power(np.subtract(prob, anal), 2))
   BS_ref = np.nanmean(np.power(np.subtract(climo, anal), 2)) # BS of climatological forecast vs obs (climo forecast = 0, i.e, no anomaly)
   BSS = 1 - (BS/BS_ref)
   return BSS

def get_ABN(ensemble_data, clim_upper, clim_lower):
   total = np.array([[len(ensemble_data[:,r,c]) for c, column in enumerate(row)] for r, row in enumerate(ensemble_data[0])])
   above = np.divide(np.sum(ensemble_data>clim_upper, axis=0), total)
   below = np.divide(np.sum(ensemble_data<clim_lower, axis=0), total)
   normal = np.divide(np.sum((ensemble_data<clim_lower)+(ensemble_data>clim_upper)==0, axis=0), total)
   nan_mask_upper = np.isnan(clim_upper)
   nan_mask_lower = np.isnan(clim_lower)
   nan_mask_normal = np.isnan(clim_lower)+np.isnan(clim_upper)
   above[nan_mask_upper]=np.nan
   below[nan_mask_lower]=np.nan
   normal[nan_mask_normal]=np.nan
   return (above, below, normal)

def plot_temp_probability(num, model, lons, lats, data, fstart, fend, BSS_under=None, BSS_over=None, BSS_normal=None, 
                          analysis=False, init=None, domain=None, parallels=None, merid=None, latrange=None, lonrange=None, figsize=None):
   if analysis:
      print(f'Plot {num} ... Temperature Anomaly Probability for {str(model).upper()} (Analysis)')
   else:
      print(f'Plot {num} ... Temperature Anomaly Probability for {str(model).upper()}')
   print('========================================')
   print('Prepping Basemap')
   # Get the basemap object
   if str(domain).upper() == 'OTHER':
      proj, fig = get_domains(domain=domain, parallels=parallels, merid=merid, latrange=latrange, lonrange=lonrange, figsize=figsize, num=num)
   else:
      proj, fig = get_domains(domain=domain, num=num)
      lonrange, latrange = get_domains(domain=domain, return_box=True)
   
   ax = fig.add_axes([.05,.05,.9,.9], projection=proj)
   ax.set_extent([lonrange[0], lonrange[1], latrange[0], latrange[1]], crs=ccrs.PlateCarree())
   countries = cf.NaturalEarthFeature(
      category='cultural',
      name='admin_0_boundary_lines_land',
      scale='50m',
      facecolor='none')
   states = cf.NaturalEarthFeature(
      category='cultural',
      name='admin_1_states_provinces_lines',
      scale='50m',
      facecolor='none')
   coastlines = cf.NaturalEarthFeature(
      category='physical',
      name='coastline',
      scale='50m',
      facecolor='none')
   #ocean = cf.NaturalEarthFeature(
      #category='physical',
      #name='ocean',
      #scale='50m',
      #facecolor='none')
   land = cf.NaturalEarthFeature(
      category='physical',
      name='land',
      scale='50m',
      facecolor='none')
   lakes = cf.NaturalEarthFeature(
      category='physical',
      name='lakes',
      scale='50m',
      facecolor='none')
   ax.add_feature(land, zorder=2)
   ax.add_feature(cf.OCEAN, zorder=100, edgecolor='white', facecolor='white')
   ax.add_feature(lakes, edgecolor='black', facecolor='none', linewidth=.5, zorder=2)
   ax.add_feature(countries, edgecolor='black', facecolor='none', linewidth=.7, zorder=4)
   ax.add_feature(states, edgecolor='black', facecolor='none', linewidth=.5, zorder=6)
   ax.add_feature(coastlines, edgecolor='black', facecolor='none', linewidth=.4, zorder=5)
   latlongrid=5.
   parallels = list(np.arange(0., 90., latlongrid))
   meridians = list(np.arange(180., 360., latlongrid))
   
   print('Basemap Prepped')

   print('Plotting Temperature Error ...')
   clevs = [-1., -.9, -.8, -.7, -.6, -.5, -.4, -.3333, .3333, .4, .5, .6, .7, .8, .9, 1.]
   n=len(clevs)-1
   vmin, vmax = (-1., 1.)
   MEG_error_colors = [
      '#1b2c62', # -1. - -.9 
      '#214c93', # -.9 - -.8 
      '#3b7bbd', # -.8 - -.7
      '#66afde', # -.7 - -.6 
      '#94d3f3', # -.6 - -.5
      '#c1e7f8', # -.5 - -.4 
      '#e6f5fb', # -.4 - -.33
      '#A9A9A9', # -.33 - .33 #ffffff is white and #D3D3D3 is very light gray
      '#fdf3c5', # .33 - 0.4 
      '#fdd778', # 0.4 - 0.5 
      '#fdaa31', # 0.5 - 0.6 
      '#f77c2b', # 0.6 - 0.7 
      '#e44529', # 0.7 - 0.8 
      '#c31d24', # 0.8 - 0.9 
      '#921519', # 0.9 - 1.0

   ]
   cmap = colors.ListedColormap(MEG_error_colors)
   norm = colors.BoundaryNorm(clevs, cmap.N, clip=False)

   #plot_temp_prob = ax.pcolormesh(lons, lats, data, shading='flat', cmap=cmap, vmin=vmin, vmax=vmax, zorder=1, transform=ccrs.PlateCarree(), norm=norm)
   #Use contourf when using finer_grid
   plot_temp_prob = ax.contourf(lons, lats, data, shading='flat', levels=clevs, cmap=cmap, vmin=vmin, vmax=vmax, zorder=1, transform=ccrs.PlateCarree(), norm=norm)
   print('Plotted')

   valid_start = init+td(hours=fstart)
   valid_end = init+td(hours=fend)

   model_name = model.upper().replace('_', ' ')
   if valid_start.year == valid_end.year:
      if valid_start.month == valid_end.month:
         valid_string = valid_start.strftime(' valid %HZ %d - ') + valid_end.strftime('%d %B %Y')
      else:
         valid_string = valid_start.strftime(' valid %HZ %d %B - ') + valid_end.strftime('%d %B %Y')
   else:
      valid_string = valid_start.strftime(' valid %HZ %d %B %Y - ') + valid_end.strftime('%HZ %d %B %Y')
   if analysis:
      title_string_left = f'{model_name} Analysis |' + valid_string
   else:
      title_string_left = f'{model_name} | ' + init.strftime('initialized %HZ %d %B %Y |') + valid_string
   title_string_right = r'2-m T Anomaly Probability (%)' 
   plt.title(title_string_left, loc='left', fontsize=7)
   plt.title(title_string_right, loc='right', fontsize=7)

   ax = plt.gca()
   if not analysis:
      above_text = ax.text(x=1.11, y=0.85, s='Above\nNormal\nBSS: {:.2f}'.format(BSS_over), transform=ax.transAxes, horizontalalignment='center', verticalalignment='center', color='black', fontsize=10, zorder=10) 
      below_text = ax.text(x=1.11, y=0.15, s='Below\nNormal\nBSS: {:.2f}'.format(BSS_under), transform=ax.transAxes, horizontalalignment='center', verticalalignment='center', color='black', fontsize=10, zorder=10) 
      normal_text = ax.text(x=1.11, y=0.5, s='Near\nNormal\nBSS: {:.2f}'.format(BSS_normal), transform=ax.transAxes, horizontalalignment='center', verticalalignment='center', color='black', fontsize=10, zorder=10) 
   else:
      above_text = ax.text(x=1.11, y=0.85, s='Above\nNormal', transform=ax.transAxes, horizontalalignment='center', verticalalignment='center', color='black', fontsize=10, zorder=10) 
      below_text = ax.text(x=1.11, y=0.15, s='Below\nNormal', transform=ax.transAxes, horizontalalignment='center', verticalalignment='center', color='black', fontsize=10, zorder=10) 
      normal_text = ax.text(x=1.11, y=0.5, s='Near\nNormal', transform=ax.transAxes, horizontalalignment='center', verticalalignment='center', color='black', fontsize=10, zorder=10) 

   fig.subplots_adjust(left=.05, right=.91, top=.95, bottom=.05, wspace=0, hspace=0)
   divider = make_axes_locatable(ax)
   cbax = divider.new_horizontal(size="3%", pad=0.1, axes_class=plt.Axes)
   cax = fig.add_axes(cbax)
   #cax = fig.add_axes([.91, .05, .01, .9])
   cbar_ticks = [-1., -.9, -.8, -.7, -.6, -.5, -.4, -.3333, .3333, .4, .5, .6, .7, .8, .9, 1.]
   cb = plt.colorbar(plot_temp_prob, orientation='vertical', cax=cax, cmap=cmap, norm=norm, boundaries=clevs, spacing='uniform', ticks=cbar_ticks, drawedges=True)
   cb.dividers.set_color('black')
   cb.dividers.set_linewidth(2)
   cb.ax.tick_params(size=0, labelsize=8, labelright=True, labelleft=False, right=False)
   cb.ax.set_yticklabels(['100','90','80','70','60','50','40','33','33','40','50','60','70','80','90','100'])
   cax.hlines([0, 1], 0, 1, colors='black', linewidth=4)
    
   if analysis:
      save_string = f'{str(model).upper()}_ANALYSIS_2MTEMP_ANOM_PROBABILITY_{domain}.'+valid_start.strftime('valid_%Y%m%d%H')+valid_end.strftime('_%Y%m%d%H.png')
   else:
      save_string = f'{str(model).upper()}_2MTEMP_ANOM_PROBABILITY_{domain}.'+init.strftime('init_%Y%m%d%H.')+valid_start.strftime('valid_%Y%m%d%H')+valid_end.strftime('_%Y%m%d%H.png')
   save_path = os.path.join(SAVE_DIR, save_string)
   plt.savefig(save_path, bbox_inches='tight')
   print(f'Plot {num} saved successfully in {save_path}')
   plt.close(num)
   print('========================================')

def main(lon_range=None, lat_range=None):

   if str(analysis).upper() == 'RAP':
      analysis_string = rap_string
   elif str(analysis).upper() == 'GFS':
      analysis_string = gefs_string # GEFS produces "gfs" member data files
   elif str(analysis).upper() == 'GEFS':
      analysis_string = gefs_string
   # Get, Process, and Plot Data
   # Each dataset needs to be processed separately because they contain different precip accumulation step ranges
   num = 0
   # Given init_date
   fstart = day_range[0]*24
   fend = day_range[1]*24+24   
   for domain in domains:
      if str(domain).upper() != 'OTHER':
         lon_range, lat_range = get_domains(domain=domain, return_box=True)
      for init in inits:
         for model in models:
               # if a file isn't found, we'll ignore this init and go to the next one
               try: 
                  if str(model).upper() == 'GEFS':
                     gefs_ens_num = 30
                     gefs_ens_ctl = 'c00'
                     gefs_ens_mems = ['p{:02d}'.format(num) for num in np.arange(1, gefs_ens_num+1)]

                     # Get analysis temperature data
                     if str(analysis).upper() == 'GEFS':
                        use_mem_for_anal = 'ctl'
                     elif str(analysis).upper() == 'GFS':
                        use_mem_for_anal = 'gfs'
                     elif str(analysis).upper() == 'RAP':
                        use_mem_for_anal = 'rap'
                     else:
                        raise OSError(f'Script does not support using {analysis} as the analysis dataset, use \'gfs\', \'rap\', or \'gefs\'')

                     anal_temp, anal_y, anal_x = ens_lower_data(init, str(analysis).upper(), analysis_string, fstart, fend, ens_mem=None, analysis=True)
                     climo_upper_temp, climo_lower_temp, climo_lats, climo_lons = climo_lower_data(init, climo_string, fstart, fend)

                     if str(analysis).upper() in ['GFS','GEFS']:
                        anal_lons, anal_lats = np.meshgrid(anal_x, anal_y)
                     climo_lons, climo_lats = np.meshgrid(climo_lons, climo_lats)

                     if finer_grid:
                        if (anal_y[0] > anal_y[-1]).all():
                           lats_fine = np.arange(anal_y[0], anal_y[-1], -1*finer_grid)
                        else:
                           lats_fine = np.concatenate((np.arange((anal_y[0]).all(), np.min(anal_y), -1*finer_grid), np.arange(np.max(anal_y), (anal_y[-1]).all(), -1*finer_grid)))
                        if (anal_x[0] < anal_x[-1]).all():
                           lons_fine = np.arange(anal_x[0], anal_x[-1], finer_grid)
                        else:
                           lons_fine = np.concatenate((np.arange((anal_x[0]).all(), np.max(anal_x), finer_grid), np.arange(np.min(anal_x), (anal_x[-1]).all(), finer_grid)))
                        lons_fine, lats_fine = np.meshgrid(lons_fine, lats_fine)
                        anal_temp = diff_data(lons_fine, lats_fine, anal_x, anal_y, anal_temp)
                        use_lons, use_lats = (lons_fine, lats_fine)
                     else:
                        use_lons, use_lats = (anal_x, anal_y)

                     gefs_temps = []
                     for p, gefs_ens_mem in enumerate(gefs_ens_mems):
                        try:
                           gefs_mem_temp, gefs_y, gefs_x = ens_lower_data(init, str(model).upper(), gefs_string, fstart, fend, ens_mem=gefs_ens_mem)
                           if p == 0:
                              gefs_lons, gefs_lats = np.meshgrid(gefs_x, gefs_y)
                           gefs_mem_temp_interp = diff_data(use_lons, use_lats, gefs_lons, gefs_lats, gefs_mem_temp)
                           gefs_temps.append(gefs_mem_temp_interp)
                        except OSError:
                           continue

                     gefs_temps = np.array(gefs_temps)
                     climo_upper_temp_interp = diff_data(use_lons, use_lats, climo_lons, climo_lats, climo_upper_temp)
                     climo_lower_temp_interp = diff_data(use_lons, use_lats, climo_lons, climo_lats, climo_lower_temp)

                     gefs_over, gefs_under, gefs_normal = get_ABN(gefs_temps, climo_upper_temp_interp, climo_lower_temp_interp)
                     anal_over, anal_under, anal_normal = get_ABN(np.array([anal_temp]), climo_upper_temp_interp, climo_lower_temp_interp)
                     
                     gefs_heatmap = np.multiply(np.nanmax([gefs_under, gefs_normal, gefs_over], axis=0), np.argmax([gefs_under, gefs_normal, gefs_over], axis=0)-1)
                     anal_heatmap = np.multiply(np.nanmax([anal_under, anal_normal, anal_over], axis=0), np.argmax([anal_under, anal_normal, anal_over], axis=0)-1) 
                     #Original code
                     #use_lons, use_lats = (gefs_lons, gefs_lats)

                     #mask1 = gefs_lons < lon_range[0]
                     #mask2 = gefs_lons > lon_range[1]
                     #mask3 = gefs_lats > lat_range[1]
                     #mask4 = gefs_lats < lat_range[0]

                     #Set if using finer_grid
                     mask1 = lons_fine < lon_range[0]
                     mask2 = lons_fine > lon_range[1]
                     mask3 = lats_fine > lat_range[1]
                     mask4 = lats_fine < lat_range[0]

                     #Set if finer_grid=None (using anal grid)
                     #mask1 = anal_x < lon_range[0]
                     #mask2 = anal_x > lon_range[1]
                     #mask3 = anal_y > lat_range[1]
                     #mask4 = anal_y < lat_range[0]
                     mask = mask1+mask2+mask3+mask4
                     mgefs_over, mgefs_under, mgefs_normal = [np.ma.MaskedArray(gefs_prob, mask=mask) for gefs_prob in [gefs_over, gefs_under, gefs_normal]]
                     manal_over, manal_under, manal_normal = [np.ma.MaskedArray(anal_prob, mask=mask) for anal_prob in [anal_over, anal_under, anal_normal]]
                     
                     BSS_under = BSS(mgefs_under, manal_under)               
                     BSS_over = BSS(mgefs_over, manal_over)               
                     BSS_normal = BSS(mgefs_normal, manal_normal)               

                     if str(domain).upper() == 'OTHER':
                        plot_temp_probability(num, model, use_lons[1:-1,1:-1], use_lats[1:-1,1:-1], gefs_heatmap[1:-1,1:-1], fstart, fend, BSS_under=BSS_under, BSS_over=BSS_over, BSS_normal=BSS_normal, init=init, domain=domain, parallels=parallels, merid=merid, lonrange=lon_range, latrange=lat_range, figsize=figsize)
                     else:
                        plot_temp_probability(num, model, use_lons[1:-1,1:-1], use_lats[1:-1,1:-1], gefs_heatmap[1:-1,1:-1], fstart, fend, BSS_under=BSS_under, BSS_over=BSS_over, BSS_normal=BSS_normal, init=init, domain=domain)
                     num+=1
                     if str(domain).upper() == 'OTHER':
                        plot_temp_probability(num, analysis, use_lons[1:-1,1:-1], use_lats[1:-1,1:-1], anal_heatmap[1:-1,1:-1], fstart, fend, analysis=True, init=init, domain=domain, parallels=parallels, merid=merid, lonrange=lon_range, latrange=lat_range, figsize=figsize)
                     else:
                        plot_temp_probability(num, analysis, use_lons[1:-1,1:-1], use_lats[1:-1,1:-1], anal_heatmap[1:-1,1:-1], fstart, fend, analysis=True, init=init, domain=domain)
                     num+=1
               except OSError as e:
                  print(e)
                  print("Continuing ...")
                  continue

main()
