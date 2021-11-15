# ==========================================
# Title: Height Anomalies
# Description: Plot anomaly statistics using grib output of
# or several models, obs, and climatology (GP height contours)
# Author: Marcel Caron
# Created: Aug. 31, 2021
# Last Modified: Oct. 27, 2021
# ==========================================

import os
import numpy as np
import cfgrib
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cf
from datetime import datetime, timedelta as td
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
from scipy.constants import knot
import seaborn as sns
import pickle

# ============= CHANGES BELOW ==============

DATA_DIR = '/scratch2/NCEPDEV/ovp/Marcel.Caron/MEG/data/'
SAVE_DIR = '/scratch2/NCEPDEV/ovp/Marcel.Caron/MEG/cases/tests/'

# ~~~~~~~
# This is the continuation of DATA_DIR, because each model directory is unique.
# Use standard C datetime format codes (e.g., %Y for year) for init datetimes, 
# and use VALIDY and VALIDM for valid year and month, respectively
# ~~~~~~~
climo_string = 'noscrub/climo_files/mean_VALIDMVALIDD' 
cfsanl_string = 's2s/cfs.VALIDYVALIDM/pgbhnl.gdas.VALIDYVALIDM.grib2'
cfs_string = 's2s/cfs.%Y%m%d/pgbf.01.%Y%m%d%H.VALIDYVALIDM.avrg.grib.%HZ.grb2' 
gefs_string = 'common/gefs.%Y%m%d/geavg.t%Hz.pgrb2a.0p50.fFFF'
gefsanl_string = 'common/gefs.%Y%m%d/gegfs.t%Hz.pgrb2a.0p50.f000'
gfsanl_string = 'common/gfs.%Y%m%d/gfs.t%Hz.pgrb2.0p25.f00'

# ~~~~~~~
# Valid times
# Can be a single datetime or a tuple of datetimes. In the latter case, plots 
# will be a mean value of the period requested
# ~~~~~~~
valids = [(datetime(2021,1,16,12,0,0), datetime(2021,1,21,12,0,0)),
          (datetime(2021,1,18,12,0,0), datetime(2021,1,25,12,0,0))] # inclusive

# ~~~~~~~
# step between each valid in the period requested (in hours)
# must include this if any valids are tuples, otherwise this is ignored
# ~~~~~~~
step = 12 

# ~~~~~~~
# e.g., ['gefs', 'cfs']
# ~~~~~~~
models = ['gefs'] 

# ~~~~~~~
# Analysis field used for comparison
# Options are gefs, cfs, gfs
# ~~~~~~~
analysis = 'gefs'

# ~~~~~~~
# Create a plots for each model at the following init time
# ~~~~~~~
init = datetime(2021,1,10,12,0,0)

# ~~~~~~~
# e.g., [1000, 700, 500, 250]
# ~~~~~~~
plevels=[500,700]

# ~~~~~~~
# list of domains to plot
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

# tuple width, height for desired dimensions of plots
figsize = (9.5,7.) 

# ~~~~~~~
# Defined requested grid resolution in degrees, or put None
# ~~~~~~~
finer_grid = None

# ============= CHANGES ABOVE ==============


# input types are datetime, datetime, timedelta objects
def daterange(start, end, td):
   curr = start
   while curr <= end:
      yield curr
      curr+=td

# Works generally if you're plotting in the western hemisphere
def get_lons(da):
   if np.max(da.longitude) > 180.:
      lons = np.subtract(
         np.mod(np.add(np.array(da.longitude),180.),360.),180.
      )
   else:
      lons = np.array(da.longitude)
   return lons

# Give a string, output is either the map projection and figure objects or 
# lonrange and latrange tuples that define the desired domain
# String can specify a predefined domain, in which case settings are defined in 
# the function below, or create a user-defined domain if the string == 'other'
# Feel free to edit the settings if needed (increasing figsize width can 
# adjust placement of colorbar, for example
def get_domains(domain='northamer', parallels=None, merid=None, latrange=None, 
                lonrange=None, num=0, figsize=None, return_box=False):
   if str(domain).upper() == 'NORTHAMER':
      if not return_box:
         fig = plt.figure(num=num, figsize=(7.6,7.))
      parallels = (11., 31.)
      merid = -100.
      latrange = (8., 72.)
      lonrange = (-151., -49.)
   if str(domain).upper() == 'FRONTRANGE':
      if not return_box:
         fig = plt.figure(num=num, figsize=(10.,7.))
      parallels = (30., 35.)
      merid = -97.
      latrange = (36., 46.)
      lonrange = (-112.5, -96.5)
   elif str(domain).upper() == 'NWHEAT':
      if not return_box:
         fig = plt.figure(num=num, figsize=(7.7,7.))
      parallels = (25., 35.)
      merid = -122.5
      latrange = (42., 49.1)
      lonrange = (-127.3, -116.3)
   elif str(domain).upper() == 'CONUS':
      if not return_box:
         fig = plt.figure(num=num, figsize=(9.,7.))
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
   proj = ccrs.LambertConformal(
      central_longitude=merid, standard_parallels=parallels
   )
   if return_box:
      return lonrange, latrange
   else:
      return proj, fig

def get_leadhour(init, valid):
   return str(int((valid-init).days*24. + (valid-init).seconds/3600.))

def upper_data(init, valid, model, model_string, level, step=None, 
               use_valid=False):
   
   print('========================================')
   print(f'Getting and Summing {model} Geopotential Height Data')
   
   if isinstance(valid, tuple):
      if use_valid:
         fpaths = [
            os.path.join(
               DATA_DIR, 
               val.strftime(model_string)
                  .replace('VALIDY', val.strftime('%Y'))
                  .replace('VALIDM', val.strftime('%m'))
                  .replace('VALIDD', val.strftime('%d'))
                  .replace('FFF', get_leadhour(init, val))
            ) 
            for val in daterange(valid[0], valid[1], td(hours=step))
         ]
      else:
         fpaths = [
            os.path.join(
               DATA_DIR, 
               init.strftime(model_string)
                  .replace('VALIDY', val.strftime('%Y'))
                  .replace('VALIDM', val.strftime('%m'))
                  .replace('VALIDD', val.strftime('%d'))
                  .replace('FFF', get_leadhour(init, val))
            ) 
            for val in daterange(valid[0], valid[1], td(hours=step))
         ]
   else:
      if use_valid:
         fpaths = [
            os.path.join(
               DATA_DIR, 
               valid.strftime(model_string)
                  .replace('VALIDY', valid.strftime('%Y'))
                  .replace('VALIDM', valid.strftime('%m'))
                  .replace('VALIDD', valid.strftime('%d'))
                  .replace('FFF', get_leadhour(init, valid))
            )
         ]
      else:
         fpaths = [
            os.path.join(
               DATA_DIR, 
               init.strftime(model_string)
                  .replace('VALIDY', valid.strftime('%Y'))
                  .replace('VALIDM', valid.strftime('%m'))
                  .replace('VALIDD', valid.strftime('%d'))
                  .replace('FFF', get_leadhour(init, valid))
            )
         ]
   climo_indices = {'1000':0, 
                    '700':1,
                    '500':2,
                    '250':3}
   cfs_indices = {'1000':0,
                    '975':1,
                    '950':2,
                    '925':3,
                    '900':4,
                    '875':5,
                    '850':6,
                    '825':7,
                    '800':8,
                    '775':9,
                    '750':10,
                    '700':11,
                    '650':12,
                    '600':13,
                    '550':14,
                    '500':15,
                    '450':16,
                    '400':17,
                    '350':18,
                    '300':19,
                    '250':20,
                    '225':21,
                    '200':22,
                    '175':23}

   gefs_indices = {'1000':0,
                    '925':1,
                    '850':2,
                    '700':3,
                    '500':4,
                    '250':6,
                    '200':7}
   gfs_indices = {'1000':0,
                    '975':1,
                    '950':2,
                    '925':3,
                    '900':4,
                    '850':6,
                    '800':8,
                    '750':10,
                    '700':11,
                    '650':12,
                    '600':13,
                    '550':14,
                    '500':15,
                    '450':16,
                    '400':17,
                    '350':18,
                    '300':19,
                    '250':20,
                    '200':22}

   if str(model).upper()[:4] == 'GEFS':
      model_indices = gefs_indices
   elif str(model).upper()[:3] == 'CFS':
      model_indices = cfs_indices
   elif str(model).upper()[:3] == 'GFS':
      model_indices = gfs_indices
   elif str(model).upper() == 'CLIMO':
      pass
   else:
      print(f'Warning: Model name doesn\'t match models for which this code is'
            + f' configured; consider developing new settings for the'
            + f' requested model: {str(model).upper()}')
   if model_string[-3:] == '.nc':
      data = []
      for f, fpath in enumerate(fpaths):
         with Dataset(fpath, mode='r') as fh:
            if f==0:
               lons = fh.variables['lon'][:]
               lats = fh.variables['lat'][:]
            data.append(fh.variables['HGT'][:])
      data = np.nanmean(data, axis=0)
   else:
      data = []
      for f, fpath in enumerate(fpaths):
         type_of_level = 'isobaricInhPa'
         short_name = 'gh'
         isobaric_indexing=True
         with xr.open_dataset(fpath, engine='cfgrib', backend_kwargs=dict(
                  indexpath='',
                  filter_by_keys={'typeOfLevel':type_of_level,
                                  'shortName':short_name}
                  )) as ds_gh:
               if str(model).upper() == 'CLIMO':
                  if isobaric_indexing:
                     # in the narrow use case I'm interested in, climo has a 
                     # both time and vertical dimension
                     hgt = ds_gh.gh[climo_indices[str(level)]] 
                  else:
                     # in the narrow use case I'm interested in, climo has a 
                     # both time and vertical dimension
                     hgt = ds_gh.gh 
               else:
                  if isobaric_indexing:
                     # may need to be changed if adapting this script to use 
                     # other model inputs
                     hgt = ds_gh.gh[model_indices[str(level)]] 
                  else:
                     hgt = ds_gh.gh
               if f==0:
                  lons = get_lons(hgt)
                  lats = np.array(hgt.latitude)
               data.append(np.array(hgt))
      data = np.nanmean(data, axis=0)
   
   # convert gpm to dam
   data/=10. 
   
   print('Data Retrieved ' + u'\u2713')
   print('========================================') 
   return (data, lats, lons)

# Returns numpy data meshgrid as a bilinear interpolation of model data onto 
# analysis grid
def diff_data(anl_lons, anl_lats, model_lons, model_lats, model_data):
   anl_lons, anl_lats, model_lons, model_lats, model_data = [
      np.array(item) 
      for item in [anl_lons, anl_lats, model_lons, model_lats, model_data]
   ]
   model_data_interpd = griddata(
      (model_lons.flatten(), model_lats.flatten()),
      model_data.flatten(),
      (anl_lons, anl_lats), method='linear')
   return model_data_interpd

def plot_anom(num, model, lons, lats, data, valid, plevel, init=None, 
              domain=None, lonrange=None, latrange=None, parallels=None, 
              merid=None):
   print(f'Plot {num} ... {plevel} hPa Heights Anomaly plot for '
         + f'{str(model).upper()}')
   print('========================================')
   print('Prepping Basemap')
   if str(domain).upper() == 'OTHER':
      proj, fig = get_domains(
         domain=domain, parallels=parallels, merid=merid, latrange=latrange, 
         lonrange=lonrange, figsize=figsize, num=num
      )
   else:
      proj, fig = get_domains(domain=domain, num=num)
      lonrange, latrange = get_domains(domain=domain, return_box=True) 
   
   ax = fig.add_axes([.05,.05,.9,.9], projection=proj)
   ax.set_extent(
      [lonrange[0], lonrange[1], latrange[0], latrange[1]], 
      crs=ccrs.PlateCarree()
   )
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
   ocean = cf.NaturalEarthFeature(
      category='physical',
      name='ocean',
      scale='50m',
      facecolor='none')
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
   ax.add_feature(ocean, zorder=2)
   ax.add_feature(
      lakes, edgecolor='black', facecolor='none', linewidth=.5, zorder=2
   )
   ax.add_feature(
      countries, edgecolor='black', facecolor='none', linewidth=.7, zorder=4
   )
   ax.add_feature(
      states, edgecolor='black', facecolor='none', linewidth=.5, zorder=6
   )
   ax.add_feature(
      coastlines, edgecolor='black', facecolor='none', linewidth=.4, zorder=5
   )
   latlongrid=5.
   parallels = list(np.arange(0., 90., latlongrid))
   meridians = list(np.arange(180., 360., latlongrid))
   
   print('Basemap Prepped ' + u'\u2713')
   print('Plotting Height Anomalies ...')
   clevs = [-999., -14., -12., -10., -8., -6., -4., -2., 
            2., 4., 6., 8., 10., 12., 14., 999.]
   MEG_anom_colors = [
      '#1b2c62', # < -14. F
      '#214c93', # -14. - -12. F 
      '#3b7bbd', # -12. - -10. F 
      '#66afde', # -10. - -8.0 F 
      '#94d3f3', # -8.0 - -6.0 F 
      '#c1e7f8', # -6.0 - -4.0 F 
      '#e6f5fb', # -4.0 - -2.0 F 
      '#ffffff', # -2.0 - 2.00 F 
      '#fdf3c5', # 2.00 - 4.00 F 
      '#fdd778', # 4.00 - 6.00 F 
      '#fdaa31', # 6.00 - 8.00 F 
      '#f77c2b', # 8.00 - 10.0 F 
      '#e44529', # 10.0 - 12.0 F 
      '#c31d24', # 12.0 - 14.0 F 
      '#921519', # > 14.0 F  
   ]
   cmap = colors.ListedColormap(MEG_anom_colors)
   norm = colors.BoundaryNorm(clevs, cmap.N, clip=False)

   plot_anom = ax.pcolormesh(
      lons, lats, data, shading='flat', cmap=cmap, zorder=1, 
      transform=ccrs.PlateCarree(), norm=norm,
   )
  
   print('Plotted ' + u'\u2713')
   model_name = model.upper().replace('_', ' ')
   if init:
      if isinstance(valid, tuple):
         title_string_left = (f'{model_name} | ' 
                              + init.strftime('initialized %HZ %d %B %Y') 
                              + valid[0].strftime(' valid %HZ %d %B %Y') 
                              + valid[1].strftime(' to %HZ %d %B %Y'))
      else:
         title_string_left = (f'{model_name} | ' 
                              + init.strftime('initialized %HZ %d %B %Y') 
                              + valid.strftime(' valid %HZ %d %B %Y'))
   else:
      if isinstance(valid, tuple):
         title_string_left = (f'{model_name} | ' 
                              + valid[0].strftime(' valid %HZ %d %B %Y') 
                              + valid[1].strftime(' to %HZ %d %B %Y'))
      else:
         title_string_left = (f'{model_name} | ' 
                              + valid.strftime(' valid %HZ %d %B %Y'))
   title_string_right = f'{plevel}-hPa g (gpm)'
   plt.title(title_string_left, loc='left', fontsize=7)
   plt.title(title_string_right, loc='right', fontsize=7)

   fig.subplots_adjust(
      left=.05, right=.93, top=.95, bottom=.05, wspace=0.05, hspace=0
   )
   cax = fig.add_axes([.93, .05, .01, .9])
   cbar_ticks = clevs[1:-1]
   cb = plt.colorbar(
      plot_anom, orientation='vertical', cax=cax, cmap=cmap, norm=norm, 
      boundaries=clevs, spacing='uniform', ticks=cbar_ticks, drawedges=True
   )
   cb.dividers.set_color('black')
   cb.dividers.set_linewidth(2)
   cb.ax.tick_params(
      labelsize=8, labelright=True, labelleft=False, right=False
   )
   cb.ax.set_yticklabels(['-14','-12','-10','-8','-6','-4','-2',
                          '2','4','6','8','10','12','14',''])
   cax.hlines([0, 1], 0, 1, colors='black', linewidth=4)

   if init:
      if isinstance(valid, tuple):
         save_string = (f'{model}_{str(plevel).upper()}HGT_ANOMALY_{domain}.'
                        + init.strftime('init_%Y%m%d%H.')
                        + valid[0].strftime('valid_%Y%m%d%H')
                        + valid[1].strftime('_%Y%m%d%H.png'))
      else:
         save_string = (f'{model}_{str(plevel).upper()}HGT_ANOMALY_{domain}.'
                        + init.strftime('init_%Y%m%d%H.')
                        + valid.strftime('valid_%Y%m%d%H.png'))
   else:
      if isinstance(valid, tuple):
         save_string = (f'{model}_{str(plevel).upper()}HGT_ANOMALY_{domain}.'
                        + valid[0].strftime('valid_%Y%m%d%H')
                        + valid[1].strftime('_%Y%m%d%H.png'))
      else:
         save_string = (f'{model}_{str(plevel).upper()}HGT_ANOMALY_{domain}.'
                        + valid.strftime('valid_%Y%m%d%H.png'))
   save_path = os.path.join(SAVE_DIR, save_string)
   plt.savefig(save_path, bbox_inches='tight')
   print(f'Plot {num} saved successfully in {save_path}')
   plt.close(num)
   print('========================================')

def plot_height(num, model, lons, lats, data_hgt, valid, plevel, 
                contour_every=6, domain=None, parallels=None, merid=None, 
                latrange=None, lonrange=None, figsize=None, init=None):
   print(f'Plot {num} ... {plevel} hPa Height plot for {str(model).upper()}')
   print('========================================')
   print('Prepping Basemap')
   # Get the basemap object
   if str(domain).upper() == 'OTHER':
      proj, fig = get_domains(
         domain=domain, parallels=parallels, merid=merid, latrange=latrange, 
         lonrange=lonrange, figsize=figsize, num=num
      )
   else:
      proj, fig = get_domains(domain=domain, num=num)
      lonrange, latrange = get_domains(domain=domain, return_box=True) 

   ax = fig.add_axes([.05,.05,.9,.9], projection=proj)
   ax.set_extent(
      [lonrange[0], lonrange[1], latrange[0], latrange[1]], 
      crs=ccrs.PlateCarree()
   )
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
   ocean = cf.NaturalEarthFeature(
      category='physical',
      name='ocean',
      scale='50m',
      facecolor='#5c5c5c')
   land = cf.NaturalEarthFeature(
      category='physical',
      name='land',
      scale='50m',
      facecolor='#d3d3d3')
   lakes = cf.NaturalEarthFeature(
      category='physical',
      name='lakes',
      scale='50m',
      facecolor='#5c5c5c')
   ax.add_feature(land, zorder=2)
   ax.add_feature(ocean, zorder=2)
   ax.add_feature(lakes, edgecolor='black', linewidth=.5, zorder=2)
   ax.add_feature(
      countries, edgecolor='black', facecolor='none', linewidth=.7, zorder=4
   )
   ax.add_feature(
      states, edgecolor='black', facecolor='none', linewidth=.5, zorder=6
   )
   ax.add_feature(
      coastlines, edgecolor='black', facecolor='none', linewidth=.4, zorder=5
   )
   latlongrid=5.
   parallels = list(np.arange(0., 90., latlongrid))
   meridians = list(np.arange(180., 360., latlongrid))

   print('Basemap Prepped ' + u'\u2713')

   # really rough approximations using hypsometric equation
   contour_low = np.round(675.18*np.log(1000./float(plevel))) 
   contour_high = np.round(867.06*np.log(1000./float(plevel)))

   if contour_low%2==1:
      # no odd contours please
      contour_low+=1 
   contour_ints = np.arange(contour_low,contour_high,contour_every)
   # Mask height data that are out of bounds; helps clabel better-populate the 
   # domain with contour labels
   # If the contours don't fill the plot adequately, it might help to increase 
   # the buffer (+/- deg. lon/lat) on either side of the mask
   mask1 = lons < lonrange[0] - 15. 
   mask2 = lons > lonrange[1] + 15.
   mask3 = lats > latrange[1] - 15.
   mask4 = lats < latrange[0] + 15.
   mask = mask1+mask2+mask3+mask4 
   mdata_hgt = np.ma.MaskedArray(data_hgt,mask=mask)

   print('Plotting Heights ...')
   plot_hgt = ax.contour(
      lons, lats, mdata_hgt, levels=contour_ints, colors='black', 
      linewidths=1.5, transform=ccrs.PlateCarree(), zorder=9
   )
   plt.clabel(
      plot_hgt, contour_ints, fmt='%i', fontsize=8.5, inline=1, 
      inline_spacing=0,  manual=False
   )

   print('Plotted ' + u'\u2713')
   model_name = model.upper().replace('_', ' ')
   if init:
      if isinstance(valid, tuple):
         title_string_left = (f'{model_name} | ' 
                              + init.strftime('initialized %HZ %d %B %Y') 
                              + valid[0].strftime(' valid %HZ %d %B %Y') 
                              + valid[1].strftime(' to %HZ %d %B %Y'))
      else:
         title_string_left = (f'{model_name} | ' 
                              + init.strftime('initialized %HZ %d %B %Y') 
                              + valid.strftime(' valid %HZ %d %B %Y'))
   else:
      if isinstance(valid, tuple):
         title_string_left = (f'{model_name} | ' 
                              + valid[0].strftime(' valid %HZ %d %B %Y') 
                              + valid[1].strftime(' to %HZ %d %B %Y'))
      else:
         title_string_left = (f'{model_name} | ' 
                              + valid.strftime(' valid %HZ %d %B %Y'))
   title_string_right = f'{plevel}-hPa g (dam)'
   plt.title(title_string_left, loc='left', fontsize=7)
   plt.title(title_string_right, loc='right', fontsize=7)

   if init:
      if isinstance(valid, tuple):
         save_string = (f'{model}_Z{str(plevel).upper()}_{domain}.'
                        + init.strftime('init_%Y%m%d%H.')
                        + valid[0].strftime('valid_%Y%m%d%H')
                        + valid[1].strftime('_%Y%m%d%H.png'))
      else:
         save_string = (f'{model}_Z{str(plevel).upper()}_{domain}.'
                        + init.strftime('init_%Y%m%d%H.')
                        + valid.strftime('valid_%Y%m%d%H.png'))
   else:
      if isinstance(valid, tuple):
         save_string = (f'{model}_Z{str(plevel).upper()}_{domain}.'
                        + valid[0].strftime('valid_%Y%m%d%H')
                        + valid[1].strftime('_%Y%m%d%H.png'))
      else:
         save_string = (f'{model}_Z{str(plevel).upper()}_{domain}.'
                        + valid.strftime('valid_%Y%m%d%H.png'))
   save_path = os.path.join(SAVE_DIR, save_string)
   plt.savefig(save_path, bbox_inches='tight')
   print(f'Plot {num} saved successfully in {save_path}')
   plt.close(num)
   print('========================================')

# Get, Process, and Plot Data
def main(lonrange=None, latrange=None):

   num = 0
   if str(analysis).upper() == 'CFS':
      analysis_string = cfsanl_string
      use_valid = False
   elif str(analysis).upper() == 'GEFS':
      analysis_string = gefsanl_string
      use_valid = True
   elif str(analysis).upper() == 'GFS':
      analysis_string = gfsanl_string
      use_valid = True
   for domain in domains:
      if str(domain).upper() != 'OTHER':
         # these variables don't exist unless we have defined the "other" 
         # domain
         lonrange, latrange = get_domains(domain=domain, return_box=True) 
      for model in models:
         if str(model).upper() == 'CFS':
            fcst_string = cfs_string
         elif str(model).upper() == 'GEFS':
            fcst_string = gefs_string
         for valid in valids:
            for plevel in plevels:
               # if a file isn't found, we'll ignore this init and go to the 
               # next one 
               try: 
                  fcst_hgt, fcst_lats, fcst_lons = upper_data(
                     init, valid, str(model).upper(), fcst_string, plevel, 
                     step=step
                  )
                  obs_hgt, obs_y, obs_x = upper_data(
                     init, valid,  str(analysis).upper()+' ANALYSIS', 
                     analysis_string, plevel, step=step, use_valid=use_valid
                  )
                  climo_hgt, climo_lats, climo_lons = upper_data(
                     init, valid, 'CLIMO', climo_string, plevel, 
                     use_valid=use_valid, step=step
                  )
                  fcst_lons, fcst_lats = np.meshgrid(fcst_lons, fcst_lats)
                  obs_lons, obs_lats = np.meshgrid(obs_x, obs_y)
                  climo_lons, climo_lats = np.meshgrid(climo_lons, climo_lats)
                  climo_interp = diff_data(
                     obs_lons, obs_lats, climo_lons, climo_lats, climo_hgt
                  )
                  fcst_interp = diff_data(
                     obs_lons, obs_lats, fcst_lons, fcst_lats, fcst_hgt
                  )
                  obs_anom = np.subtract(obs_hgt, climo_interp)
                  fcst_anom = np.subtract(fcst_interp, climo_interp)
                  anom_lons, anom_lats = (obs_lons, obs_lats)
                  if finer_grid:
                     if obs_y[0] > obs_y[-1]:
                        lats_fine = np.arange(
                           obs_y[0], obs_y[-1], -1*finer_grid
                        )
                     else:
                        lats_fine = np.concatenate((
                           np.arange(obs_y[0], np.min(obs_y), -1*finer_grid), 
                           np.arange(np.max(obs_y), obs_y[-1], -1*finer_grid)
                        ))
                     if obs_x[0] < obs_x[-1]:
                        lons_fine = np.arange(obs_x[0], obs_x[-1], finer_grid)
                     else:
                        lons_fine = np.concatenate((
                           np.arange(obs_x[0], np.max(obs_x), finer_grid), 
                           np.arange(np.min(obs_x), obs_x[-1], finer_grid)
                        ))
                     lons_fine, lats_fine = np.meshgrid(lons_fine, lats_fine)
                     fcst_anom = diff_data(
                        lons_fine, lats_fine, obs_lons, obs_lats, fcst_anom
                     )
                     obs_anom = diff_data(
                        lons_fine, lats_fine, obs_lons, obs_lats, obs_anom
                     )
                     anom_lons, anom_lats = (lons_fine, lats_fine)
                  if str(domain).upper() == 'OTHER':
                     plot_height(
                        num, str(model).upper(), fcst_lons[1:-1,1:-1], 
                        fcst_lats[1:-1,1:-1], fcst_hgt[1:-1,1:-1], valid, 
                        plevel, init=init, domain=domain, lonrange=lonrange, 
                        latrange=latrange, parallels=parallels, merid=merid
                     )
                  else:
                     plot_height(
                        num, str(model).upper(), fcst_lons[1:-1,1:-1], 
                        fcst_lats[1:-1,1:-1], fcst_hgt[1:-1,1:-1], valid, 
                        plevel, init=init, domain=domain
                     )
                  num+=1
                  if str(domain).upper() == 'OTHER':
                     plot_height(
                        num, str(analysis).upper()+' ANALYSIS', 
                        obs_lons[1:-1,1:-1], obs_lats[1:-1,1:-1], 
                        obs_hgt[1:-1,1:-1], valid, plevel, init=None, 
                        domain=domain, lonrange=lonrange, latrange=latrange, 
                        parallels=parallels, merid=merid
                     )
                  else:
                     plot_height(
                        num, str(analysis).upper()+' ANALYSIS', 
                        obs_lons[1:-1,1:-1], obs_lats[1:-1,1:-1], 
                        obs_hgt[1:-1,1:-1], valid, plevel, init=None, 
                        domain=domain
                     )
                  num+=1
                  if str(domain).upper() == 'OTHER':
                     plot_height(
                        num, 'CLIMO', climo_lons[1:-1,1:-1], 
                        climo_lats[1:-1,1:-1], climo_hgt[1:-1,1:-1], valid, 
                        plevel, init=None, domain=domain, lonrange=lonrange, 
                        latrange=latrange, parallels=parallels, merid=merid
                     )
                  else:
                     plot_height(
                        num, 'CLIMO', climo_lons[1:-1,1:-1], 
                        climo_lats[1:-1,1:-1], climo_hgt[1:-1,1:-1], valid, 
                        plevel, init=None, domain=domain
                     )
                  num+=1
                  if str(domain).upper() == 'OTHER':
                     plot_anom(
                        num, str(model).upper(), anom_lons[1:-1,1:-1], 
                        anom_lats[1:-1,1:-1], fcst_anom[1:-1,1:-1], valid, 
                        plevel, init=init, domain=domain, lonrange=lonrange, 
                        latrange=latrange, parallels=parallels, merid=merid
                     )
                  else:
                     plot_anom(
                        num, str(model).upper(), anom_lons[1:-1,1:-1], 
                        anom_lats[1:-1,1:-1], fcst_anom[1:-1,1:-1], valid, 
                        plevel, init=init, domain=domain
                     )
                  num+=1
                  if str(domain).upper() == 'OTHER':
                     plot_anom(
                        num, str(analysis).upper(), anom_lons[1:-1,1:-1], 
                        anom_lats[1:-1,1:-1], obs_anom[1:-1,1:-1], valid, 
                        plevel, init=None, domain=domain, lonrange=lonrange, 
                        latrange=latrange, parallels=parallels, merid=merid
                     )
                  else:
                     plot_anom(
                        num, str(analysis).upper(), anom_lons[1:-1,1:-1], 
                        anom_lats[1:-1,1:-1], obs_anom[1:-1,1:-1], valid, 
                        plevel, init=None, domain=domain
                     )
                  num+=1
               except OSError as e:
                  print(e)
                  print("Continuing ...")
                  continue

main()
