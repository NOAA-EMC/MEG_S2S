# S2S\_Python
##### Contributors:  Marcel Caron \(NOAA/EMC\),
##### Description:   Python library containing tools for case-by-case 
#####                seasonal and subseasonal verification in the EMC Model 
#####                Evaluation Group.
##### Created:       Sept. 29, 2021
##### Last Modified: Oct. 25, 2021

### SET UP
---

Scripts in this repository are designed to run on NCEP's Hera supercomputer.
Via Hera, it is suggested to have the following lines in your .bashrc or .cshrc file:

```
> module use -a /contrib/hpss/modulefiles
> module load hpss/4.1.0.3
> module use -a /contrib/intel/modulefiles
> module load intel/18.0.5.274
> module use -a /contrib/anaconda/modulefiles
> module load anaconda/latest
> module use -a /contrib/met/modulefiles
> module load met/9.1
> module use -a /contrib/METplus/modulefiles
> module load metplus/3.1
> module load nco/4.9.1
```

If you use Bash:
```
> export PATH="/scratch2/NCEPDEV/ovp/Marcel.Caron/anaconda3/bin:$PATH"
> export PROJ_LIB="/scratch2/NCEPDEV/ovp/Marcel.Caron/anaconda3/share/proj"
```
If you use Csh:
```
> setenv PATH="/scratch2/NCEPDEV/ovp/Marcel.Caron/anaconda3/bin:$PATH"
> setenv PROJ_LIB="/scratch2/NCEPDEV/ovp/Marcel.Caron/anaconda3/share/proj"
```

Then on the command line:
```
$ source ~/.bashrc
```
OR
```
$ source ~/.cshrc
```

Finally, the cartopy module uses natural earth files to draw polylines.  Some python scripts in this repository use cartopy and require these natural earth files to plot features like country and continent borders.
To use these files, copy the cartopy folder and its contents from the repository directory into your home directory, like this:

```
$ mkdir -p ~/.local/share && cp -r cartopy ~/.local/share/.
```

---
### DATA DIRECTORIES

Scripts use monthly cfs and analysis timeseries, but some scripts may be flexible enough to use hourly "time" timeseries
Both monthly and hourly CFS and CFSR input data, and other relevant data, may be found in tarballs on HPSS \(as of Sept. 29 2021\):

CFS:
```
/NCEPPROD/hpssprod/runhistory/cfsYYYY/YYYYmm/YYYYmmdd/monthly/cfs.pgbf.YYYYmmddHH.mPP.monthly.tar
/NCEPPROD/hpssprod/runhistory/cfsYYYY/YYYYmm/YYYYmmdd/time/cfs.YYYYmmddHH.mPP.time.tar 
```

CFSR:
```
/NCEPPROD/hpssprod/runhistory/cfsYYYY/YYYYmm/monthly_analysis/cfs.pgb.YYYYmm.monthly.tar
/NCEPPROD/hpssprod/runhistory/cfsYYYY/YYYYmm/YYYYmmdd/Analysis/cfs.pgrbh.YYYYmmdd.tar
```

GEFS:
```
/NCEPPROD/2year/hpssprod/runhistory/rhYYYY/YYYYmm/YYYYmmdd/com_gefs_prod_gefs.YYYYmmdd_HH.atmos_pgrb2ap5.tar
/NCEPPROD/2year/hpssprod/runhistory/rhYYYY/YYYYmm/YYYYmmdd/com_gefs_prod_gefs.YYYYmmdd_HH.atmos_pgrb2sp25.tar
```

GFS:
```
/NCEPPROD/hpssprod/runhistory/rhYYYY/YYYYmm/YYYYmmdd/com_gfs_prod_gfs.YYYYmmdd_HH.gfs_pgrb2.tar
```

RAP:
```
/NCEPPROD/hpssprod/runhistory/2year/rhYYYY/YYYYmm/YYYYmmdd/com_rap_prod_rap.YYYYmmdd\(00-05/06-11/12-17/18-23\).awp130.tar
```

`YYYY` = year  
`mm`   = month  
`dd`   = day  
`HH`   = hour  
`PP`   = model member in {01, 02, 03, 04}  
For RAP, `\(string1/string2/string3/string4\)` means choose and insert a string out of {string1, string2, string3, string4}  

---
### EXAMPLE FOR DOWNLOADING DATA FROM HPSS

Example in Bash:
```
$ savedir=${savedir:-/path/to/save/directory}
$ gefstarpath=${gefsdir:-/NCEPPROD/2year/hpssprod/runhistory/rh2021/202101/20210110/com_gefs_prod_gefs.20210110_12.atmos_pgrb2ap5.tar}
$ gefsfile="atmos/pgrb2ap5/geavg.t12z.pgrb2a.0p50.f138"
$ htar -xvf ${gefstarpath} ./${gefsfile}; mkdir ${savedir}/gefs.20210110; mv ${gefsfile} ${savedir}/gefs.20210110/.
```

---
### INSTRUCTIONS FOR SETTING UP METPLUS CONFIG FILES

Note: The repo contains sample configuration files for generating s2s statistics. These files are located in the METplus folder 
For CFS: MET and METplus files are preset for isobaric temperature verification.
For GEFS: MET and METplus files are present for accumulated precipitation verification.

CFS\_TMP\_grid\_stat\_metplus.conf:
- Under \[config\], set INIT\_BEG and INIT\_END to the same desired initialization date
- Under \[config\], set BOTH\_VAR1\_LEVELS to comma-separated list of pressure levels \(e.g., P850 for 850-hPa level\)
- Under \[config\], set the LEAD\_SEQ to comma-separated list of desired lead days \(number of days followed by "d"\)
- Under \[config\], set OBTYPE to the name of the observation dataset you are using
- Under \[config\], set GRID\_STAT\_CONFIG\_FILE to path to MET config file in your system \(provided in this repo\)
- Under \[dir\], set CONFIG\_DIR, FCST\_\*\_INPUT\_DIR, OBS\_\*\_INPUT\_DIR, and GRID\_STAT\_OUTPUT\_DIR to desired directories \(INPUT\_BASE and OUTPUT\_BASE are set in hera.conf and are usable variables\) 
- Under \[dir\], set GRID\_STAT\_CLIMO\_MEAN\_INPUT\_DIR to directory containing climatological data
- Under \[filename\_templates\], set FCST\_REGRID\_DATA\_PLANE\_INPUT\_TEMPLATE and OBS\_GRID\_STAT\_INPUT\_TEMPLATE, matching names of forecast and observational data files, respectively. Use METplus datetime keywords where necessary \(e.g., {init?fmt=%Y%m%d%H}\)
- Under \[filename\_templates\], set GRID\_STAT\_CLIMO\_MEAN\_INPUT\_TEMPLATE matching names of climatological data files

GEFS\_PCP\_grid\_stat\_metplus.conf:
- Under \[config\], set VALID\_BEG and VALID\_END to the same desired valid date
- Under \[config\], set BOTH\_VAR1\_LEVELS to the desired accumulation \(e.g., 'A06' for 6-h accumulation; this should correspond to settings for FCST\_PCP\_COMBINE\_OUTPUT\_ACCUMS and OBS\_PCP\_COMBINE\_OUTPUT\_ACCUMS\)
- Under \[config\], set BOTH\_VAR1\_THRESH to comma-separated list of desired precipitation thresholds \(used for contingency table statistics\)
- Under \[config\], set the LEAD\_SEQ to comma-separated list of desired lead hours 
- Under \[config\], set the OBTYPE to the name of the observation dataset you are using
- Under \[config\], set REGRID\_DATA\_PLANE\_VERIF\_GRID to a sample analysis file \(will be used as reference grid for interpolation\)
- Under \[config\], set GRID\_STAT\_CONFIG\_FILE to path to MET config file in your system \(provided in this repo\)
- Under \[config\], you may adjust the GRID\_STAT\_NEIGHBORHOOD\_WIDTH and GRID\_STAT\_NEIGHBORHOOD\_SHAPE, which affects neighborhood output statistics \(e.g., nbrcts\)
- Under \[config\], set FCST\_PCP\_COMBINE\_OUTPUT\_ACCUMS and OBS\_PCP\_COMBINE\_OUTPUT\_ACCUMS to desired output accumulation time in hours, and OBS\_PCP\_COMBINE\_INPUT\_ACCUMS to the correct input accumulation time
- Under \[dir\], set CONFIG\_DIR, FCST\_\*\_INPUT\_DIR, FCST\_\*\_OUTPUT\_DIR, OBS\_\*\_INPUT\_DIR, OBS\_\*\_OUTPUT\_DIR, and/or GRID\_STAT\_OUTPUT\_DIR \(ditto for OBS\_...\) to desired directories \(INPUT\_BASE and OUTPUT\_BASE are set in hera.conf and are usable variables\)
- Under \[dir\], set GRID\_STAT\_CLIMO\_MEAN\_INPUT\_DIR to directory containing climatological data
- Under \[filename\_templates\], set FCST\_PCP\_COMBINE\_INPUT TEMPLATE and OBS\_PCP\_COMBINE\_INPUT\_TEMPLATE, matching name of forecast and observational data files, respectively. Use METplus datetime keywords where necessary \(e.g., {init?fmt=%Y%m%D%H}\)
- Under \[filename\_templates\], set GRID\_STAT\_CLIMO\_MEAN\_INPUT\_TEMPLATE matching names of climatological data files

---
### OPTIONAL INSTRUCTIONS FOR SETTING UP THE MET CONFIG FILES

- OPTIONAL: In the 'mask' dictionary, set 'poly' to a comma-separated list of strings containing paths to polylines, 
            if desired. May also leave blank.  Polylines are defined in .poly files, with the name of the polygon/domain 
            at the top \(this will be the "VX\_MASK" in stat output files\), followed by a column of lat/lon pairs 
            indicating polygon vertices. The first and last coordinates in the column should be the same \(closed polygon\).

---
### INSTRUCTIONS FOR RUNNING METPLUS

- Use earlier sections \(Instructions for setting up the METplus config files, and optional instructions for setting up the met config files\) to adjust settings according to your use case \(use, e.g., CFS_*_metplus.conf and GridStatConfig_CFS_TMP for CFS temperature verifcation\)
- Set up hera.conf INPUT\_DIR and OUTPUT\_DIR \(locations of input data and desired output directory of statistics files, respectively\)
- Set up metplus\_job.sh accordingly:
  - Set the --output option to where SBATCH output should be stored for each job. Include %j or $$ in the name so that each job is uniquely identifiable.
  - Set other SBATCH options as needed to run jobs from your account \(e.g., --account=ovp will not work if you are not in the ovp group\)
  - Set the conf\_base, which is the directory that contains the MET and METplus config files you'd like to use.
  - Set the metplus\_conf\_fname, which is the name of the metplus configuration file you'd like to run.
- To run, enter the following on the command line:
- `$ sbatch metplus_job.sh` 
- To check the status of your job, here are two things you can do:
  - `$ squeue -u First.Last #` to see if your job is still running
  - check the output file of your job \(name and location was defined in metplus\_job.sh as the --output option in the SBATCH header\) to see command-line output so far

---
### INSTRUCTIONS FOR SETTING UP PYTHON SCRIPTS

ensemble\_temperature\_probability.py:
- Set DATA\_DIR and SAVE\_DIR as high-level data directory \(to be continued in the file string setting below\) and the directory where plots should be saved, respectively.
- Set gefs\_string and climo\_string as continuations of DATA\_DIR, which in combination form the path to the gefs data files and climatology data files, respectively.
- Set values in the day\_range tuple, which represent start and end day of averaging period 
- Set models, a comma-separated list of strings representing requested models \(currently only "gefs" is supported\)
- Set analysis string, which will be used to plot a reference probability map \(currently only "gefs" is supported\)
- Set list of requested initialization times in datetime format
- Set list of domains over which to plot the map \(if predefined domains don't cover desired area, you may set the domain manually by including "other" in this list\)
- Set latrange, lonrange, parallels, merid, and figsize according to desired domain.  These will be used if "other" was included in the list of requested domains
- Manually set the grid resolution, in degrees lat/lon, or put None to use the analysis grid. Setting this may cause the script to fail if the requested resolution is too fine for the domain.

height\_anoms.py:
- Set DATA\_DIR and SAVE\_DIR as high-level data directory \(to be continued in the file string setting below\) and the directory where plots should be saved, respectively.
- Set climo\_string, gefs\_string and gefsanl\_string or cfs\_string and cfsanl\_string as continuations of DATA\_DIR, which in combination form the path to the climatology data files and forecast/analysis data files, respectively.
- Set list of requested valid times in datetime format. May alternatively put either a list of datetime tuples or a list of mixed datetime tuples and single datetime objects. For tuples, the script will compute a mean value across a valid period, where the tuple represents the start and end times \(inclusive\) of the valid period. 
- If datetime tuples are included in the valids list, set step in hours. Step is the time interval to step through the valid period. The valid field at each step will be used to compute the overall mean across the valid period.
- Set models, a comma-separated list of strings representing requested models
- Set analysis string, which will be used to plot a reference probability map \(climatology data is often the limiting factor here\).
- Set requested initialization time in datetime format
- Set list requested isobars. Plots will be created at each level if the field is available.
- Set list of domains over which to plot the map \(if predefined domains don't cover desired area, you may set the domain manually by including "other" in this list\)
- Set latrange, lonrange, parallels, merid, and figsize according to desired domain.  These will be used if "other" was included in the list of requested domains
- Manually set the grid resolution, in degrees lat/lon, or put None to use the analysis grid. Setting this may cause the script to fail if the requested resolution is too fine for the domain.

temp\_anoms.py:
- Set DATA\_DIR and SAVE\_DIR as high-level data directory \(to be continued in the file string setting below\) and the directory where plots should be saved, respectively.
- Set climo\_string, gefs\_string and gefsanl\_string or cfs\_string and cfsanl\_string as continuations of DATA\_DIR, which in combination form the path to the climatology data files and forecast/analysis data files, respectively.
- Set list of requested valid times in datetime format. May alternatively put either a list of datetime tuples or a list of mixed datetime tuples and single datetime objects. For tuples, the script will compute a mean value across a valid period, where the tuple represents the start and end times \(inclusive\) of the valid period. 
- If datetime tuples are included in the valids list, set step in hours. Step is the time interval to step through the valid period. The valid field at each step will be used to compute the overall mean across the valid period.
- Set models, a comma-separated list of strings representing requested models
- Set analysis string, which will be used to plot a reference probability map \(climatology data is often the limiting factor here\).
- Set requested initialization time in datetime format
- Set list requested isobars. Plots will be created at each level if the field is available.
- Set list of domains over which to plot the map \(if predefined domains don't cover desired area, you may set the domain manually by including "other" in this list\)
- Set latrange, lonrange, parallels, merid, and figsize according to desired domain.  These will be used if "other" was included in the list of requested domains
- Manually set the grid resolution, in degrees lat/lon, or put None to use the analysis grid. Setting this may cause the script to fail if the requested resolution is too fine for the domain.

modified\_taylor\_diagram.py:
- Set list of requested initialization times in datetime format
- Set STATS\_DIR and SAVE\_DIR as high-level directory containing metplus output \(to be continued in the file string setting below\) and the directory where plots should be saved, respectively.
- Set fstring string as continuation of STATS\_DIR, which in combination form the path to the metplus output text files \("cnt" or "sal1l2" suffixes\).
- Set models, a comma-separated list of strings representing requested models
- Set analysis string, which will be used to find the correct stats file\(s\)
- Set values in the day\_range tuple, which represent start and end day of averaging period 
- Set step length between each valid time in the period defined in "day\_range". Used to find the correct stats files 
- Set list requested isobars. Plots will be created at each level if the field is available.
- Set requested variable to plot. Used to find the correct stats file\(s\)
- Set list of domains/polylines used to compute statistics. Used to index stats in the stats file\(s\).  Make sure this matches one of the values under the "VX\_MASK" column in the stats file\(s\).
- Set list of colors and shapes from which the script will grab to plot each requested model.  The order will correspond to the order in which models are listed.

teleconnections.py:
- Set DATA\_DIR and SAVE\_DIR as high-level data directory \(to be continued in the file string setting below\) and the directory where plots should be saved, respectively.
- Set fcst\_string, anl\_string, and climo\_string as continuations of DATA\_DIR, which in combination form the path to the forecast, analysis, and climatology data files, respectively
- Set list of requested initialization times in datetime format
- Set values in the lead\_range tuple, which represent start and end lead months for which indices will be calculated \(1 month increments\)
- Set values in the record record\_range tuple, which represent start and end datetimes of the complete record of monthly analyses that you have available
- Set list of requested teleconnections
- Set requested model \(currently only 'cfs' is tested, but others should work in theory if the data are already in monthly format\)
- Set list of desired line colors for analysis, 'model\_name', climatology forecast, and persistence forecast indices/errors \(in that order\)

---
### INSTRUCTIONS FOR RUNNING PYTHON SCRIPTS

python\_job.sh:
- Most of the above scripts can't be run on the login nodes without getting your job killed and receiving an email.
- In python\_job.sh, set the --output option to where SBATCH output should be stored for each job. Include %j or $$ in the name so that each job is uniquely identifiable.
- Set other SBATCH options as needed to run jobs from your account \(e.g., `--account=ovp` will not work if you are not in the ovp group\)
- Set the script\_dir, the directory containing the python script you'd like to run.
- Set the name of the python script you'd like to run \(part of the srun command\)
- To run any of the above python scripts, adjust the settings in python\_job.sh as described above, then enter the following on the command line:
- `$ sbatch python_job.sh` 

