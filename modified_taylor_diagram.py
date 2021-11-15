########################################################
#
# Developed: May 3, 2021 by Marcel Caron 
# Last Modified: Oct. 27, 2021 by Marcel Caron             
# Title: Taylor Diagram, modified
# Descr.: Plots Normalized Standardized Deviation of 
# the Forecast Anomaly as a function of Anomaly 
# Correlation Coefficient, with contoured Normalized 
# Root Mean Squared Anomaly Error, in polar coordinates
#
########################################################

import matplotlib.pyplot as plt
import matplotlib.markers as mmarkers
from scipy.interpolate import interp1d
import numpy as np
import pandas as pd
from datetime import datetime, timedelta as td
import os

# ==================CHANGE BELOW=========================

# ~~~~~~~
# Create plots for each model, plev, and domain at each of the following inits
# ~~~~~~~
inits = [datetime(2021,1,10,12,0,0)]

# ~~~~~~~
# Directory containing GridStat continuous statistics ('cnt') and scalar L1L2 
# statistics ('sal1l2')
# ~~~~~~~
STATS_DIR = '/scratch2/NCEPDEV/stmp1/Shannon.Shields/scripts/s2s/cases/met_out/GEFS/'

# ~~~~~~~
# Desired directory for saving output plots
# ~~~~~~~
SAVE_DIR = '/scratch2/NCEPDEV/ovp/Marcel.Caron/MEG/cases/tests/'

# ~~~~~~~
# Continuation of STATS_DIR, in order to get the paths to statistics files
# use {model} as model name, {varbl} as variable name, {analysis} as analysis 
# name, {lead} as lead hour, {type} as statistics type ...
# standard C datetime format codes (e.g., %Y for year) for init datetimes
# and VALIDY, VALIDm, VALIDd, and VALIDH for valid year, month, day, and hour 
# ~~~~~~~
fstring = '%Y%m/grid_stat/grid_stat_{model}_{varbl}_vs_{analysis}_ANL_{varbl}_P850_{lead}0000L_VALIDYVALIDmVALIDd_VALIDH0000V_{type}.txt'

# ~~~~~~~
# Create plots for each init, plev, and domain with each of the following 
# models (will be capitalized automatically)
# ~~~~~~~
models = ['gefs'] 

# ~~~~~~~
# Name of the analysis used to compute statistics (will be capitalized 
# automatically)
# ~~~~~~~
analysis = 'gefs'

# ~~~~~~~
# Lead day range to plot, following CPC standards (will be used to compute 
# valid times, so make sure it matches stat file metadata)
# ~~~~~~~
day_range = (6, 10)

# ~~~~~~~
# step length between each valid time in the period defined in "day_range", 
# in hours
# ~~~~~~~
step = 24 

# ~~~~~~~
# levels to plot (e.g., ['2m', '500', '700'])
# ~~~~~~~
plevs = ['850']

# ~~~~~~~
# variable to plot (as defined in your METplus statistics output)
# ~~~~~~~
varbl = 'TMP'

# ~~~~~~~ 
# Create plots for each init, plev, and model at each of the following domains
# Needs to be defined in the METplus statistics output under 'VXMASK' column
# ~~~~~~~
domains = ['conus']

# ~~~~~~~
# list of colors and shapes used to plot each requested model onto the diagram.
# The lists can be any length > len(models)
# Colors and shapes will be assigned to models based on the orders listed here
# Use any color/shape format recognizable by matplotlib modules
# ~~~~~~~
colors = ['cyan', 'black', 'mediumorchid', 'red', 'gold', 'lime'] 
shapes = ['o', 'P', 'v', 's', 'd','^', '*', 'X', 'p'] 

# ==================CHANGE ABOVE=========================

def mscatter(plot_object, x, mlist):
   if (mlist is not None) and (len(mlist)==len(x)):
      paths = []
      for marker in mlist:
         if isinstance(marker, mmarkers.MarkerStyle):
            marker_obj = marker
         else:
            marker_obj = mmarkers.MarkerStyle(marker)
         path = marker_obj.get_path().transformed(marker_obj.get_transform())
         paths.append(path)
      plot_object.set_paths(paths)
      return plot_object

def plot_datapoint(ds_l1l2, df_cnt, colors_list, markers_list):
   strd_anom = np.sqrt(df_l1l2['FFABAR'][0]/df_l1l2['OOABAR'][0])   
   plt.errorbar(
      np.arccos(df_cnt['ANOM_CORR'][0]), strd_anom, 
      xerr=[
         [
            np.arccos(df_cnt['ANOM_CORR_BCU'][0])
            - np.arccos(df_cnt['ANOM_CORR'][0])
         ], 
         [
            np.arccos(df_cnt['ANOM_CORR'][0])
            - np.arccos(df_cnt['ANOM_CORR_BCL'][0])
         ]
      ], 
      color=colors_list[0], marker=markers_list[0], markersize=9, mec='black', 
      mew=0.8, ecolor='black', elinewidth=1.6, ls='none'
   )

def taylor_diagram(df_cnt, df_l1l2, save_path, figsize=(9, 8), 
                   xlabel='Anomaly Correlation Coefficient', 
                   ylabel='Normalized Standard Deviation', 
                   x_ticks=np.arange(0, 2.1, 0.1), 
                   y_ticks=np.arange(0, 1.1, 0.1), dpi=300, csi_cmap='Blues', 
                   rmse_label='Normalized Root Mean Square Error', 
                   title_left='Taylor Diagram', 
                   title_right='Unknown Threshold', num=0):

   fig, ax = plt.subplots(
      figsize=figsize, num=num, subplot_kw=dict(projection='polar')
   )

   curve = [[0, .5*np.pi],[1., 1.]]
   curvex = np.linspace(curve[0][0], curve[0][1], 500)
   curvey = interp1d(curve[0], curve[1])(curvex)
   arc_line = plt.plot(curvex, curvey, color='k', ls='solid', lw=1.)
   obs_point = plt.scatter([0], [1.], color='k', marker='.', s=100.)

   grid_ticks_r = np.arange(0, 2.01, 0.001)
   grid_ticks_theta = np.arccos(np.arange(0, 1.001, 0.001))
   r_g, theta_g = np.meshgrid(grid_ticks_r, grid_ticks_theta)
   nrmse = np.sqrt(
      np.subtract(
         np.add(np.power(r_g, 2), 1), 
         2*np.multiply(r_g, np.cos(theta_g))
      )
   )
   nrmse_contour = plt.contour(
      theta_g, r_g, nrmse, 
      [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3.], 
      colors='gray', linestyles='dashed', linewidths=1.
   )
   nrmse_fill = plt.contourf(
      theta_g, r_g, nrmse, [0, 1.0], colors=['lightgray','white'], zorder=0
   )
   plt.clabel(
      nrmse_contour, fmt='%1.2f', inline_spacing=0., 
      manual=[
         (np.arccos(0.995), 1.2), 
         (np.arccos(0.997), 1.5), 
         (np.arccos(0.9975), 1.75), 
         (np.arccos(0.95), 1.9), 
         (np.arccos(0.8), 1.9), 
         (np.arccos(0.6), 1.9), 
         (np.arccos(0.37), 1.9), 
         (np.arccos(0.12), 1.9)
      ]
   )

   df_cnt_mean = df_cnt.mean(axis=0)
   df_l1l2_mean = df_l1l2.mean(axis=0)
   
   # BCMSE is the centered MSE according to METplus statisticians
   FA_variance = np.subtract(
      df_l1l2_mean['FFABAR'], np.power(df_l1l2_mean['FABAR'],2)
   )
   OA_variance = np.subtract(
      df_l1l2_mean['OOABAR'], np.power(df_l1l2_mean['OABAR'],2)
   )
   nrmse = np.sqrt(np.divide(df_cnt_mean['BCMSE'], OA_variance))

   # Calculating the Standardized Anomaly Bias Error (i.e., Mean Difference in 
   # Standardized Forecast and Observed Anomalies)
   cbar = np.subtract(df_cnt_mean['OBAR'], df_l1l2_mean['OABAR']) 
   cvar = np.subtract(
             np.add(
                np.power(
                   np.add(
                      df_cnt_mean['OBAR'], 
                      df_l1l2_mean['OABAR']
                      ), 
                   2
                   ), 
                np.power(df_cnt_mean['OSTDEV'], 2)
                ), 
             np.power(cbar, 2)
             )

   std_fa = np.divide(df_l1l2_mean['FABAR'], np.sqrt(cvar))  
   std_oa = np.divide(df_l1l2_mean['OABAR'], np.sqrt(cvar)) 
   std_err = np.subtract(std_fa, std_oa) 

   df_cnt_mean['NRMSE'] = nrmse
   nstdev = (np.array(np.sqrt(np.divide(FA_variance, OA_variance))))
   df_cnt_mean['NSTDEV'] = nstdev
   # values tend to be ~.1
   df_cnt_mean['STDERR'] = np.divide(std_err, .1) 

   nstdev_i = df_cnt_mean['NSTDEV']
   plt.errorbar(
      np.arccos(df_cnt_mean['ANOM_CORR']), nstdev_i, 
      xerr=[
         [
            np.arccos(df_cnt_mean['ANOM_CORR'])
            - np.arccos(df_cnt_mean['ANOM_CORR_BCU'])
         ], 
         [
            np.arccos(df_cnt_mean['ANOM_CORR_BCL'])
            - np.arccos(df_cnt_mean['ANOM_CORR'])
         ]
      ], 
      color='black', marker='o', markersize=9, mec='black', mew=0.8, 
      ecolor='black', elinewidth=1.6, ls='none'
   )

   ax.set_thetamin(0)
   ax.set_thetamax(90)

   trans, _, _ = ax.get_xaxis_text1_transform(-6)
   ax.text(
      np.pi*.25, -.1, xlabel, transform=trans, 
      rotation=np.rad2deg(np.pi*.25)-90., fontsize=12, ha='center', 
      va='center'
   )

   # The "xlabel" function labels the radial axis, so use our ylabel variable 
   # as inputs here
   plt.xlabel(ylabel, fontsize=12, labelpad=25.) 
   plt.ylabel(ylabel, fontsize=12, labelpad=30.)
   theta_labels = [
      '{:.2f}'.format(theta_label) 
      for theta_label in np.cos(x_ticks)
   ]
   r_labels = np.array(['{:.1f}'.format(y_tick) for y_tick in y_ticks])
   plt.xticks(x_ticks, theta_labels)
   plt.yticks(y_ticks, r_labels)
   ax.tick_params(
      labelleft=True, labelright=True, labelbottom=True, labeltop=False
   )

   plt.title(title_left, fontsize=9, fontweight='bold', loc='left', pad=30.)
   plt.title(title_right, fontsize=9, fontweight='bold', loc='right', pad=30.)

   f = lambda m, c, ls, mew: plt.plot(
      [], [], marker=m, mec='k', mew=mew, color=c, ls=ls, lw=2.
   )[0]

   plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
   print(u'\u2713' + f'Plot {num} saved successfully in {save_path}')
   plt.close(num)

def main():
   num = 0
   if str(varbl).upper() in ['TMP', 'TEMP', 'T']:
      varbl_title = 'T'
   elif str(varbl).upper() in ['HGT', 'HEIGHT', 'Z']:
      varbl_title = 'g'
   else:
      print(f'Script hasn\'t been updated to accommodate the requested'
         + f' variable: {str(varbl).upper()}')
      print('Code edit is recommended.')
      print('Quitting ...')
      return None
   for init in inits:
      print('========================================')
      print('Gathering Data ...')
      df_cnt_set = False
      df_l1l2_set = False
      for m, model in enumerate(models):
         print(f'... for model #{m} ({str(model).upper()}) ...')
         model = str(model).upper().replace('V','v')
         fstart = day_range[0]*24
         fend = day_range[1]*24 + 24
         for l, lead in enumerate(np.arange(fstart, fend, step)):
            try:
               valid_date = init+td(hours=int(lead))
               fname_cnt = init.strftime(
                  fstring
                     .replace('VALIDY', valid_date.strftime('%Y'))
                     .replace('VALIDm', valid_date.strftime('%m'))
                     .replace('VALIDd', valid_date.strftime('%d'))
                     .replace('VALIDH', valid_date.strftime('%H'))
                     .replace('{model}', str(model))
                     .replace('{varbl}', str(varbl).upper())
                     .replace('{type}', str('cnt'))
                     .replace('{analysis}', str(analysis).upper())
                     .replace('{lead}', str(lead))
               )
               fname_l1l2 = valid_date.strftime(
                  fstring
                     .replace('VALIDY', valid_date.strftime('%Y'))
                     .replace('VALIDm', valid_date.strftime('%m'))
                     .replace('VALIDd', valid_date.strftime('%d'))
                     .replace('VALIDH', valid_date.strftime('%H'))
                     .replace('{model}', str(model))
                     .replace('{varbl}', str(varbl).upper())
                     .replace('{type}', str('sal1l2'))
                     .replace('{analysis}', str(analysis).upper())
                     .replace('{lead}', str(lead))
               )
               fpath_cnt = os.path.join(STATS_DIR, fname_cnt)
               fpath_l1l2 = os.path.join(STATS_DIR, fname_l1l2)
               df_cnt_temp = pd.read_csv(fpath_cnt, delim_whitespace=True)
               df_l1l2_temp = pd.read_csv(fpath_l1l2, delim_whitespace=True)
               if not df_cnt_set:
                  df_cnt = df_cnt_temp
                  df_cnt_set = True
               else:
                  df_cnt = pd.concat([df_cnt, df_cnt_temp])
               if not df_l1l2_set:
                  df_l1l2 = df_l1l2_temp
                  df_l1l2_set = True
               else:
                  df_l1l2 = pd.concat([df_l1l2, df_l1l2_temp])
            except OSError as e:
               print(e)
               print('Continuing ...')
               continue
         print(u'\u2713' + f' ... {str(model).upper()} statistics are'
               + f' gathered')

         if not df_l1l2_set or not df_cnt_set:
            raise ValueError(
               "Dataframe variables are undefined.  fstring may not be defined"
               + " correctly, or else the desired input stat files don\'t "
               + "exist."
            )
 
         for domain in domains:
            df_cnt_domain = df_cnt[
               df_cnt['VX_MASK'].str[:] == str(domain).upper()
            ]
            df_l1l2_domain = df_l1l2[
               df_l1l2['VX_MASK'].str[:] == str(domain).upper()
            ]
            for plev in plevs:
               if str(plev).upper() in ['2MT', '2M', 'T2M']:
                  plev_compare = '2'
                  plev_title = '2-m'
               else:
                  plev_compare = str(plev)
                  plev_title = f'{plev}-hPa'
                  
               df_cnt_plev = df_cnt_domain[
                  df_cnt_domain['FCST_LEV'].str[1:] == plev_compare
               ]
               df_l1l2_plev = df_l1l2_domain[
                  df_l1l2_domain['FCST_LEV'].str[1:] == plev_compare
               ]
               title_left = init.strftime(
                  f'{str(model).upper()} | {domain} | '
                  + f'{day_range[0]}-{day_range[1]} Day {plev_title} '
                  + f'{varbl_title} | init. %HZ %d %B %Y'
               )
               title_right = f''

               save_name = init.strftime(
                  f'MODIFIED_TAYLOR_DIAGRAM_{str(domain).upper()}_'
                  + f'{day_range[0]}-{day_range[1]}_DAY_'
                  + f'{str(plev_title).upper()}_{str(varbl).upper()}_'
                  + f'init_%Y%m%d%H.png'
               )
               save_path = os.path.join(SAVE_DIR, save_name)
               print(
                  init.strftime(
                     f'Creating plot {num} for {str(model).upper()} '
                     + f'{str(domain).upper()} {day_range[0]}-{day_range[1]} '
                     + f'Day {str(plev_title).upper()} {str(varbl).upper()}, '
                     + f'initialized %HZ %B %d, %Y ...'
                  )
               )
               taylor_diagram(
                  df_cnt_plev, df_l1l2_plev,save_path, figsize=(9, 8), 
                  xlabel='Anomaly Correlation Coefficient', 
                  ylabel='Normalized Standard Deviation of Anomalies', 
                  x_ticks = np.arccos([.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.99]), 
                  y_ticks = np.arange(0.2, 2., 0.2), dpi=300, 
                  rmse_label='Normalized Root Mean Square Error', 
                  title_left=title_left, title_right=title_right, num=num
               )
               num+=1

main()
