# ==========================================
# Title: Teleconnections Error
# Description: Calculate teleconnection indices using extended 
# Reanalyses, Monthly average forecasts, and climatology data,
# and plot index error as a function of forecast lead month for 
# Model forecast, climatological forecast, and persistence forecast
# Author: Marcel Caron
# Created: Aug. 31, 2021
# Last Modified: Oct. 27, 2021
# ==========================================
import os
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta as td
from dateutil.relativedelta import relativedelta as rd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import griddata
from scipy.stats import zscore
from sklearn.decomposition import PCA, FactorAnalysis
import sys
import time
import select
import pickle

# ==================CHANGE BELOW=========================

DATA_DIR = '/scratch2/NCEPDEV/'
SAVE_DIR = '/scratch2/NCEPDEV/stmp1/Marcel.Caron/cases/test'
# ~~~~~~~
# This is the continuation of DATA_DIR, because each model directory is unique
# Use standard C datetime format codes (e.g., %Y, %m for valid year, month) for 
# valid datetimes, and use {INITYm} for init yearmonth
# ~~~~~~~
fcst_string = 'stmp1/Marcel.Caron/data/NAO/cfs.{INITYm}/01/pgbf.01.{INITYm}0100.%Y%m.avrg.grib.00Z.grb2' 
anl_string = 'stmp1/Marcel.Caron/data/NAO/cfs.%Y%m/pgbhnl.gdas.%Y%m.grib2'
climo_string = 'ovp/Marcel.Caron/MEG/data/noscrub/climo_files/cmean_1d.1959%m15'

# ~~~~~~~
# Create plots at the following init times
# because this is based on monthly forecasts, it is recommended that only "1" 
# is used for the day to avoid errors related to edge cases
# ~~~~~~~
inits = [datetime(2020,10,1,0),datetime(2020,11,1,0)] 

# ~~~~~~~
# Lead time range in months, inclusive
# ~~~~~~~
lead_range = (0, 9) 

# ~~~~~~~
# What are the valid start and end dates of the monthly analyses archive?
# The record should extend at least 5 years, but at least 10 years is 
# recommended, and 20+ is ideal
# ~~~~~~~
record_range = (datetime(2011,1,1,0), datetime(2021,8,1,0)) 

# ~~~~~~~
# Create plots for each of the following teleconnections
# e.g., ['PNA', 'NAO', 'AO', 'AAO']
# ~~~~~~~
teleconnections = ['NAO','AO','AAO'] 

# ~~~~~~~
# e.g., 'cfs'
# ~~~~~~~
model_name = 'cfs'

# ~~~~~~~
# line colors associated with each forecast type/the analysis
# must be four colors, each corresponding to the following (in order): 
# analysis, 'model_name', climo, persistence
# ~~~~~~~
model_colors = ['#000000','#fb2020','#00dc00','#1e3cff'] 

# ==================CHANGE ABOVE=========================

# returns the settings (lon/lat boundaries and isobar) for the requested 
# teleconnection
def get_teleconnection_settings(teleconnection):
   teleconnection=str(teleconnection).upper()
   if teleconnection == 'NAO':
      ref_lev = 500
      ref_latrange = (20., 87.5)
      ref_lonrange = (-90., 60.)
   elif teleconnection == 'PNA':
      ref_lev = 500
      ref_latrange = (20., 87.5)
      ref_lonrange = (160., -60.)
   elif teleconnection == 'AO':
      ref_lev = 1000
      ref_latrange = (20., 87.5)
      ref_lonrange = (-180., 180.)
   elif teleconnection == 'AAO':
      ref_lev = 700
      ref_latrange = (-87.5, -20.)
      ref_lonrange = (-180., 180.)
   else:
      raise ValueError(f'Input teleconnection (\'{teleconnection}\') is not'
                       + f' supported.')
   return ref_lev, ref_latrange, ref_lonrange

# a datetime generator that takes start time, end time (inclusive) and a 
# timedelta or relativedelta object as input
def daterange(start, end, delta): 
   curr = start
   while curr <= end:
      yield curr
      curr+=delta

# Checks if longitude is in range(0., 360.) and if it is, convert to 
# range(-180., 180.)
def get_lons(da):
   if np.max(da.longitude) > 180.:
      lons = np.subtract(
         np.mod(np.add(np.array(da.longitude),180.),360.),180.
      )
   else:
      lons = np.array(da.longitude)
   return lons

# Open Analyses and create full dataset
# - get the field of interest depending on the teleconnection we're interested 
#   in (e.g., 500-, 700-, or 1000-hPa heights)
# - return the data for a single valid time and space coordinates
def get_data(fname, valid, level=500, dset='analysis', init=None, 
             return_coords=False):
   fpath = os.path.join(
      DATA_DIR, 
      valid
         .strftime(fname)
         .replace('{INITYm}',init.strftime('%Y%m'))
   )
   print(f'Opening {fpath} ...')
   with xr.open_dataset(
         fpath, engine='cfgrib', backend_kwargs=dict(
            indexpath='', 
            filter_by_keys={'typeOfLevel':'isobaricInhPa',
                            'shortName':'gh'}
         )) as ds:
      data = ds.loc[dict(isobaricInhPa=slice(level,level))].gh
   if len(data.shape) == 3:
      data = data[0]
   elif len(data.shape) == 4:
      data = data[0,0]
   elif len(data.shape) == 1:
      print(f'Warning: field component in the {dset} dataset is'
            + f' 1-dimensional.')
   else:
      pass
   if return_coords:
      lons = get_lons(data)
      lats = np.array(data.latitude)
      return np.array(data), lons, lats
   else:
      return np.array(data)

# Creates a timeseries based on a valid range
# essentially, "get_data" is looped over each valid time, and the output is 
# appended into a timeseries that is returned here
def create_dataset(fname, valid_start, valid_end, init=None, dset='analysis', 
                   level=500):
   valid_timeseries = []
   dataset = []
   return_coords = True
   count = 0
   for valid in daterange(valid_start, valid_end, rd(months=1)):
      try:
         if np.logical_not(return_coords):
            dataset.append(
               get_data(fname, valid, level=level, dset=dset, init=init)
            )
         else:
            data, lons, lats = get_data(
               fname, valid, level=level, dset=dset, init=init, 
               return_coords=return_coords
            )
            dataset.append(data)
            return_coords=False
         valid_timeseries.append(valid)
      except KeyError as e:
         print(e)
         print(f'Thrown when dealing with {dset} file valid on '
               + f'{valid.strftime("%d %B, %Y")}')
         print('Continuing ...')
         count+=1
      except EOFError as e:
         print(e)
         print(f'Thrown when dealing with {dset} file valid on '
               + f'{valid.strftime("%d %B, %Y")}')
         print('Continuing ...')
         count+=1
   print(f'{count} files were omitted due to either a KeyError or an '
         + f'EOFError, whereas {len(dataset)} were opened successfully')
   print(f'Dataset shape: {np.array(dataset).shape}')
   time.sleep(10)
   return np.array(dataset), lons, lats, np.array(valid_timeseries)

# Creates a persistence timeseries based on the analysis dataset
def create_persistence(anl_data, anl_time, valid_start, valid_end):
   time_arg_start = np.argwhere(anl_time >= valid_start)[0,0]
   time_arg_end = np.argwhere(anl_time >= valid_end)[0,0]
   persist_data = np.array(
      [anl_data[time_arg_start-1] 
      for item in range(time_arg_start-1,time_arg_end)]
   )

   # Return persistence forecast, time for which the persis. fcst. is valid
   return persist_data, anl_time[time_arg_start:time_arg_end+1] 

# Slice the analyses depending on the domain of interest (which varies by 
# teleconnection) 
def get_domain_slice(data, lats, lons, lat_range, lon_range):
   data, lons = check_lon_bounds(data, lons, lat_range, lon_range)
   lat_max = np.argwhere(lats>lat_range[1])[-1,0]+1
   lat_min = np.argwhere(lats<lat_range[0])[0,0]
   if lon_range[0] > lon_range[1]:
      lon_min = np.argwhere(lons>lon_range[0])[0,0]
      lon_max = np.argwhere(lons<lon_range[1])[-1,0]+1
   else:
      if (max(lon_range)-min(lon_range)) < (max(lons)-min(lons)):
         lon_min = np.argwhere(lons<lon_range[0])[-1,0]+1
         lon_max = np.argwhere(lons>lon_range[1])[0,0]
      else:
         lon_min = 0
         lon_max = -1
   return (data[:, lat_max:lat_min, lon_min:lon_max], 
      lats[lat_max:lat_min], 
      lons[lon_min:lon_max]) 

# gridded data may have boundaries that lie somewhere in the middle of the 
# teleconnection domain, which is no good
# if that is the case, reform the data, setting the boundaries outside or on 
# the edge of the teleconnection domain
def check_lon_bounds(data, lons, lat_range, lon_range):
   reset_bounds = True
   if lon_range[0] > lon_range[1] and lons[0] < lons[-1]:
      lons_break = np.argmin(np.abs(lons))
   elif lon_range[0] < lon_range[1] and lons[0] > lons[-1]:
      lons_break = np.argmin(lons)
   else:
      reset_bounds = False
   if reset_bounds:
      data1 = data[:, lons_break:]
      data2 = data[:, :lons_break]
      lons1 = lons[lons_break:]
      lons2 = lons[:lons_break]
      data = np.hstack((data1, data2))
      lons = np.hstack((lons1, lons2))
   return data, lons

# Returns numpy data meshgrid as a bilinear interpolation of model data onto 
# analysis grid
def diff_data(anl_lons, anl_lats, model_lons, model_lats, model_data):
   if len(anl_lons.shape) < 2 and len(anl_lats.shape) < 2:
      anl_lons, anl_lats = np.meshgrid(anl_lons, anl_lats)
   if len(model_lons.shape) < 2 and len(model_lats.shape) < 2:
      model_lons, model_lats = np.meshgrid(model_lons, model_lats)
   if (len(anl_lons.shape) >= 2 and len(anl_lats.shape) >= 2) 
         and (len(model_lons.shape) >= 2 and len(model_lats.shape) >= 2):
      anl_lons, anl_lats, model_lons, model_lats, model_data = [
         np.array(item) 
         for item in [anl_lons, anl_lats, model_lons, model_lats, model_data]
      ]
      model_data_interpd = griddata(
         (model_lons.flatten(), model_lats.flatten()),
         model_data.flatten(),
         (anl_lons, anl_lats), method='linear'
      )
      return model_data_interpd
   else:
      raise ValueError('Input latitude and longitude variables are not the '
                       + 'same shape.')

# Insert the forecast timeseries into the analysis timeseries, such that 
# analysis timeseries extends from the beginning of the record 
# (record_range[0]) up to the initialization time, and forecast timeseries 
# extends from the initialization time/initial valid time (lead month=0) to 
# the final valid time (lead month = leadrange[1]).
# Optionally, return a clipped analysis timeseries as the reference timeseries 
# (clipped where the input analysis timeseries represents a valid time that 
# is beyond the final valid time). 
def insert_subset(anl_data, subset_data, anl_time, subset_time, 
                  subset_dset='forecast', return_new_anl=True):
   if subset_time[-1] <= anl_time[-1]:
      subset_start = np.argwhere(anl_time>=subset_time[0])[0,0]
      combined = np.concatenate((anl_data[:subset_start], subset_data))
      if return_new_anl:
         new_anl_data = anl_data[:subset_start+len(subset_time)] 
         return combined, new_anl_data
      else:
         return combined
   else:
      raise ValueError(f'Subset data timeseries extends beyond the analysis '
                       + f'record (currently computing for the {subset_dset} '
                       + f'subset).')

# Return the first principal component and first eof of the principal 
# component analysis transformation of input X.
# Rotated PCA is optional (method='varimax' or method='quartimax' when running 
# the function).  
# Varimax rotated PCA is the method used most frequently because it yields 
# more physically interpretable eofs.
def get_pca(X, n_components=None, method='unrotated'):
   if str(method).upper()=='UNROTATED':
      fa = PCA()
   elif str(method).upper()=='VARIMAX':
      # will not work without scikit-learn >= v0.24.2
      fa = FactorAnalysis(rotation='varimax') 
   elif str(method).upper()=='QUARTIMAX':
      # will not work without scikit-learn >= v0.24.2
      fa = FactorAnalysis(rotation='varimax') 
   else:
      raise ValueError('Unrecognized method for pca.')
   fa.set_params(n_components=n_components)
   fa.fit(X)
   pcs = fa.transform(X)
   eofs = fa.components_
   return pcs[:,0], eofs[0]

# Get the first EOF of some dataset.
# Input data is centered, weighted by latitude, and its feature dimension 
# (lat/lon) is flattened, before it's used as input for PCA.
# Two methods for PCA are provided; the "efficient" method uses scikit-learn 
# and has the option for rotated PCA, so it is the recommended method.
def eof_analysis(data, lat, method='unrotated', efficient=True):
   data = np.array(data)
   X = data - data.mean(axis=0)
   wgts = np.sqrt(np.cos(np.deg2rad(lat)).clip(0., 1.))[..., np.newaxis]
   X_w = wgts*X
   records = X_w.shape[0]
   channels = np.product(X_w.shape[1:])
   X_f = X_w.reshape([records, channels])
   if efficient:
      pc1, E1 = get_pca(X_f, n_components=10, method=method)
   else:
      print('Warning: Efficient set to False, and as written this script can '
            + 'only compute unrotated PCA when run in inefficent mode.\n'
            + 'Continuing for unrotated PCA ...')
      A, Lh, E = np.linalg.svd(X_f, full_matrices=False)
      P = A*Lh
      pc1 = P[:,0]
      E1 = E[0]
   eof1 = np.reshape(E1, X[0].shape)
   return pc1, eof1, wgts, records, channels

# PC1 may be found for an input timeseries (e.g., forecast) by projecting the 
# dataset onto the first EOF of a reference timeseries (e.g., analysis).
# This method of determining PC1 is referred to in this script as "pobs".
def get_projected_pc1(ref_eof, data, wgts, records, channels):
   data = np.array(data)
   Y_w = wgts*data
   Y_f = Y_w.reshape([records, channels])
   pc_proj = np.dot(Y_f, ref_eof.flatten().T)
   return pc_proj

# Input data locations and other settings, output the principal components 
# determined via "pmod" and "pobs" methods (straight-up PCA and projected 
# methods, resp.), along with some time variables.
def get_all_pc1s(anl_string, fcst_string, climo_string, record_range, 
                 valid_range, init, level, lat_range, lon_range):
   # Get analysis data
   anl_data, anl_lons, anl_lats, anl_time = create_dataset(
      anl_string, record_range[0], record_range[1], init=init, dset='analysis', 
      level=level
   )
   # Ditto for cfs forecast
   fcst_data, fcst_lons, fcst_lats, fcst_time = create_dataset(
      fcst_string, valid_range[0], valid_range[1], init=init, dset='forecast', 
      level=level
   )
   # Ditto for climo forecast
   climo_data, climo_lons, climo_lats, climo_time = create_dataset(
      climo_string, valid_range[0], valid_range[1], init=init, dset='climo', 
      level=level
   )

   # Slice the spatial dimension analysis dataset based on requested domain 
   # (lat_range, lon_range)
   anl_data_sliced, anl_lats_sliced, anl_lons_sliced = get_domain_slice(
      anl_data, anl_lats, anl_lons, lat_range, lon_range
   )

   # Interpolate the forecasts to the analysis slice 
   fcst_data_interpd = np.array([
      diff_data(anl_lons_sliced, anl_lats_sliced, fcst_lons, fcst_lats, field) 
      for field in fcst_data
   ])
   climo_data_interpd = np.array([
      diff_data(
         anl_lons_sliced, anl_lats_sliced, climo_lons, climo_lats, field
      ) 
      for field in climo_data
   ])

   # Determine persistence forecast (just use preceding analysis for each time 
   # step)
   persist_data, persist_time = create_persistence(
      anl_data_sliced, anl_time, valid_range[0], valid_range[1]
   )

   # Create combined analysis and subset dataset (replace datapoints 
   # representative of forecast time period with the forecast subset)
   fcst_data_full, anl_data_full = insert_subset(
      anl_data_sliced, fcst_data_interpd, anl_time, fcst_time
   )
   climo_data_full = insert_subset(
      anl_data_sliced, climo_data_interpd, anl_time, climo_time, 
      subset_dset='climo', return_new_anl=False
   )
   persist_data_full = insert_subset(
      anl_data_sliced, persist_data, anl_time, persist_time, 
      subset_dset='persistence', return_new_anl=False
   )

   # Run PCA on combined analysis dataset
   anl_pc1, anl_eof1, wgts, records, channels = eof_analysis(
      anl_data_full, anl_lats_sliced, method='varimax'
   )

   # Use projection of first EOF of analysis PCA onto forecast datasets to get 
   # projected PC1
   fcst_pc1_pobs = get_projected_pc1(
      anl_eof1, fcst_data_full, wgts, records, channels
   )
   climo_pc1_pobs = get_projected_pc1(
      anl_eof1, climo_data_full, wgts, records, channels
   )
   persist_pc1_pobs = get_projected_pc1(
      anl_eof1, persist_data_full, wgts, records, channels
   )
   fcst_pc1_pmod, fcst_eof1, fcst_wgts, fcst_records, fcst_channels = (
      eof_analysis(fcst_data_full, anl_lats_sliced, method='varimax')
   )
   climo_pc1_pmod, climo_eof1, climo_wgts, climo_records, climo_channels = (
      eof_analysis(climo_data_full, anl_lats_sliced, method='varimax')
   )

   persist_output = (
      eof_analysis(persist_data_full, anl_lats_sliced, method='varimax')
   )
   persist_pc1_pmod = persist_output[0]
   persist_eof1 = persist_output[1]
   persist_wgts = persist_output[2]
   persist_records = persist_output[3]
   persist_channels = persist_output[4]
   valid_start_idx = -1*len(fcst_time)
   return (anl_pc1, fcst_pc1_pobs, fcst_pc1_pmod, climo_pc1_pobs, 
      climo_pc1_pmod, persist_pc1_pobs, persist_pc1_pmod, fcst_time, 
      valid_start_idx)

# Given two datetimes, returns the estimated difference in months
def diff_month(valid, init):
   return (valid.year-init.year) * 12 + valid.month - init.month

# returns the "error", or differences between fcst and anl datapoints, as 
# well as the root mean squared "error".
def get_error(anl, fcst):
   error = np.subtract(fcst, anl)
   rmse = get_rmse(error)
   return np.subtract(fcst, anl), rmse

# Given "error", returns the root mean squared "error"
def get_rmse(error):
   return np.sqrt(np.sum(np.power(error,2)))  

# Gets teleconnection indices (proxied by the first principal component of the PCA of the given timeseries) for analysis, forecast, climo, and persistence forecast timeseries ("pmod" and "pobs" for the latter 3)
# ... and collects all indices, as well as time data and data metrics, into a pandas dataframe for plotting
def get_all_teleconnections(anl_string, fcst_string, climo_string, 
                            record_range, valid_range, init, teleconnections):
   
   model_types = []
   tel_types = []
   tel_pmods = []
   tel_pmod_errorses = []
   tel_pmod_rmseses = [] 
   tel_pobses = []
   tel_pobs_errorses = []
   tel_pobs_rmseses = [] 
   valid_times = []
   lead_months = []
   for tc, teleconnection in enumerate(teleconnections):
      ref_lev, ref_latrange, ref_lonrange = get_teleconnection_settings(
         teleconnection
      )
      get_all_pc1s_output = get_all_pc1s(
         anl_string, fcst_string, climo_string, record_range, valid_range, 
         init, ref_lev, ref_latrange, ref_lonrange
      )
      tel_anl = get_all_pc1s_output[0]
      tel_fcst_pobs = get_all_pc1s_output[1]
      tel_fcst_pmod = get_all_pc1s_output[2]
      tel_climo_pobs = get_all_pc1s_output[3]
      tel_climo_pmod = get_all_pc1s_output[4]
      tel_persist_pobs = get_all_pc1s_output[5]
      tel_persist_pmod = get_all_pc1s_output[6]
      fcst_time = get_all_pc1s_output[7]
      valid_start_idx = get_all_pc1s_output[8] 
      model_type = np.hstack((
         ['ANALYSIS' for item in tel_anl[valid_start_idx:]],
         [str(model_name).upper() for item in tel_fcst_pmod[valid_start_idx:]],
         ['CLIMO' for item in tel_climo_pmod[valid_start_idx:]],
         ['PERSIST' for item in tel_persist_pmod[valid_start_idx:]]
      ))

      tel_type = np.array([teleconnection for item in model_type])

      anl_idx = zscore(tel_anl)[valid_start_idx:]
      fcst_pmod_idx = zscore(tel_fcst_pmod)[valid_start_idx:]
      climo_pmod_idx = zscore(tel_climo_pmod)[valid_start_idx:]
      persist_pmod_idx = zscore(tel_persist_pmod)[valid_start_idx:]
      fcst_pobs_idx = zscore(tel_fcst_pobs)[valid_start_idx:]
      climo_pobs_idx = zscore(tel_climo_pobs)[valid_start_idx:]
      persist_pobs_idx = zscore(tel_persist_pobs)[valid_start_idx:]
    
      anl_error, anl_rmse = get_error(anl_idx, anl_idx)
      anl_rmse_list = [anl_rmse for item in anl_error]
      fcst_pmod_error, fcst_pmod_rmse = get_error(anl_idx, fcst_pmod_idx)
      fcst_pmod_rmse_list = [fcst_pmod_rmse for item in fcst_pmod_error]
      fcst_pobs_error, fcst_pobs_rmse = get_error(anl_idx, fcst_pobs_idx)
      fcst_pobs_rmse_list = [fcst_pobs_rmse for item in fcst_pobs_error]
      climo_pmod_error, climo_pmod_rmse = get_error(anl_idx, climo_pmod_idx)
      climo_pmod_rmse_list = [climo_pmod_rmse for item in climo_pmod_error]
      climo_pobs_error, climo_pobs_rmse = get_error(anl_idx, climo_pobs_idx)
      climo_pobs_rmse_list = [climo_pobs_rmse for item in climo_pobs_error]
      persist_pmod_error, persist_pmod_rmse = get_error(
         anl_idx, persist_pmod_idx
      )
      persist_pmod_rmse_list = [
         persist_pmod_rmse for item in persist_pmod_error
      ]
      persist_pobs_error, persist_pobs_rmse = get_error(
         anl_idx, persist_pobs_idx
      )
      persist_pobs_rmse_list = [
         persist_pobs_rmse for item in persist_pobs_error
      ]

      tel_pmod = np.hstack((
         anl_idx, fcst_pmod_idx, climo_pmod_idx, persist_pmod_idx
      ))
      tel_pobs = np.hstack((
         anl_idx, fcst_pobs_idx, climo_pobs_idx, persist_pobs_idx
      ))
      tel_pmod_errors = np.hstack((
         anl_error, fcst_pmod_error, climo_pmod_error, persist_pmod_error
      ))
      tel_pobs_errors = np.hstack((
         anl_error, fcst_pobs_error, climo_pobs_error, persist_pobs_error
      ))
      tel_pmod_rmses = np.hstack((
         anl_rmse_list, fcst_pmod_rmse_list, climo_pmod_rmse_list, 
         persist_pmod_rmse_list
      ))
      tel_pobs_rmses = np.hstack((
         anl_rmse_list, fcst_pobs_rmse_list, climo_pobs_rmse_list, 
         persist_pobs_rmse_list
      ))
      valid_time = np.hstack((
         fcst_time, fcst_time, fcst_time, fcst_time
      ))
      lead_month = [diff_month(valid, init) for valid in valid_time]

      model_types = np.hstack((model_types, model_type))
      tel_types = np.hstack((tel_types, tel_type))
      tel_pmods = np.hstack((tel_pmods, tel_pmod))
      tel_pmod_errorses = np.hstack((tel_pmod_errorses, tel_pmod_errors))
      tel_pmod_rmseses = np.hstack((tel_pmod_rmseses, tel_pmod_rmses))
      tel_pobses = np.hstack((tel_pobses, tel_pobs))
      tel_pobs_errorses = np.hstack((tel_pobs_errorses, tel_pobs_errors))
      tel_pobs_rmseses = np.hstack((tel_pobs_rmseses, tel_pobs_rmses))
      valid_times = np.hstack((valid_times, valid_time))
      lead_months = np.hstack((lead_months, lead_month))
            
   df = pd.DataFrame(dict(
      MODEL=model_types, TELECONNECTION=tel_types, VALID=valid_times, 
      LEAD=lead_months, PMOD=tel_pmods, POBS=tel_pobses, 
      PMOD_ERROR=tel_pmod_errorses, PMOD_RMSE=tel_pmod_rmseses, 
      POBS_ERROR=tel_pobs_errorses, POBS_RMSE=tel_pobs_rmseses
   ))

   return df

# Selects rows in a given dataframe based on the requested teleconnection
def select_df_subset(df, teleconnection='NAO'):
   df_subset = df.loc[df['TELECONNECTION'] == str(teleconnection).upper()]   
   return df_subset

# Creates and saves plots of Mean indices and errors as functions of lead 
# month (two plots) for a given teleconnection.
def plot_stats(df_subset, inits, model_colors, teleconnection='NAO', 
               save_dir='', xlabel='Lead Months', method='pmod', 
               models=['analysis', 'forecast', 'climo', 'persist'], 
               figsize=(9,8), dpi=300):
   init_strings = [init.strftime('%b_%Y') for init in inits]
   cmap = colors.ListedColormap(model_colors)

   if str(method).upper() == 'PMOD':
      values_name = 'PMOD'
      errors_name = 'PMOD_ERROR'
      rmse_name = 'PMOD_RMSE'
   elif str(method).upper() == 'POBS':
      values_name = 'POBS'
      errors_name = 'POBS_ERROR'
      rmse_name = 'POBS_RMSE'
   
   pivot_error = pd.pivot_table(
      df_subset, index='LEAD', columns='MODEL', values=errors_name
   ).fillna(method='bfill')
   pivot_values = pd.pivot_table(
      df_subset, index='LEAD', columns='MODEL', values=values_name
   ).fillna(method='bfill')
   pivot_rmses = pd.pivot_table(
      df_subset, index='LEAD', columns='MODEL', values=rmse_name
   ).fillna(method='bfill')
   pivot_rmses_list = pivot_rmses.mean(axis=0)

   ylabels=['Error','Index']
   yticks=[np.arange(-3., 3.1, .5), np.arange(-3., 3.1, .5)]
   title_left = (f'{str(teleconnection).upper()} | Monthly Initializations '
                + f'{init_strings[0].replace("_"," ")} to '
                + f'{init_strings[-1].replace("_"," ")}')
   f_leg = lambda m, c, mew: plt.plot(
      [], [], marker=m, mec='white', mew=mew, color=c, ls='-', lw=2.
   )[0]

   for m, pivot_tbl in enumerate([pivot_error, pivot_values]):
      fig, ax = plt.subplots(num=m, figsize=figsize)
      pivot_tbl.plot(
         kind='line', colormap=cmap, legend=False, marker='o', ms=6, 
         mec='white', mew=.5, lw=2, ax=ax
      )
      plt.axhline(y=0, color='black', linestyle='--', linewidth=1, zorder=0)
      ax.set_xlabel(xlabel, fontsize=12, fontweight='bold', labelpad=15.)
      ax.set_ylabel(ylabels[m], fontsize=12, fontweight='bold', labelpad=3.)
      ax.set_xticklabels(pivot_tbl.index.astype(int))
      ax.set_yticks(yticks[m])
      ax.set_xticks(pivot_tbl.index)
      ax.set_xlim(min(pivot_tbl.index)-.1, max(pivot_tbl.index)+.1)
      ax.set_ylim(min(yticks[m]), max(yticks[m]))
      ax.tick_params(
         labelleft=True, labelright=False, labelbottom=True, labeltop=False
      )
      legend_handles = [f_leg('o', c, '.5') for c in model_colors]
      legend_labels = [
         f'{str(models[r]).upper()} ({rmse:.2f})' 
         for r, rmse in enumerate(pivot_rmses_list)
      ]

      ax.legend(
         legend_handles, legend_labels, loc='center left', fontsize=9, 
         framealpha=1, bbox_to_anchor=(1, 0.5), frameon=True, numpoints=2
      )
      fig.subplots_adjust(
         left=.05, right=.95, top=.95, bottom=.05, wspace=.05, hspace=0
      )
      title_right = f'{str(ylabels[m]).upper()} ' +r'($\sigma$)'
      ax.set_title(
         title_left, loc='left', fontsize=10, fontweight='bold', pad=15.
      )
      ax.set_title(
         title_right, loc='right', fontsize=10, fontweight='bold', pad=15.
      )
      ax.grid(
         b=True, which='major', axis='both', alpha=.5, linestyle='--', 
         linewidth=.5, zorder=0
      )
      save_name = '_'.join((
         str(teleconnection).upper(), str(method).upper(), 
         str(ylabels[m]).upper(), str(pivot_tbl.index.astype(int)[0]), 
         str(pivot_tbl.index.astype(int)[-1]), 'MONTHS', 'INIT', 
         str(init_strings[0]).upper(), str(init_strings[-1]).upper()
      ))
      save_path = os.path.join(save_dir, save_name + '.png')
      fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
      plt.close(m)

def main():
   for i, init in enumerate(inits):
      valid_range = (
         init+rd(months=lead_range[0]), 
         init+rd(months=lead_range[1])
      )
      
      # update a dataframe for each requested initialization time
      df_temp = get_all_teleconnections(
         anl_string, fcst_string, climo_string, record_range, valid_range, 
         init, teleconnections
      )
      if i == 0:
         df = df_temp
      else:
         df = pd.concat([df, df_temp])

   # Plot each requested teleconnection
   for teleconnection in teleconnections:
      df_subset = select_df_subset(df, teleconnection=teleconnection)
      plot_stats(
         df_subset, inits, model_colors, teleconnection=teleconnection, 
         save_dir=SAVE_DIR, models=['analysis', model_name, 'climo', 'persist']
      )

main()
