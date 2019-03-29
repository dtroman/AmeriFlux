# Setup to calculate Direct spectral correction method from EddyPro full
# cospectra as decribed in Polonik et al., 2019 Ag & Forest Met
# Script is publically available. If used, please include citation.
# Created by Pascal Polonik and Stephen Chan, Dec 2018

# Please be aware that all cospectra in the selected time are loaded 
# simultaneously and therefore require substantial RAM and time. 
# Cospectra are saved as .mat files when running the first execution.

# -------------------------------------------------------------#
# IMPORT PACKAGES
# -------------------------------------------------------------#
# Tested with python version 3.7
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os
import DirectCorrectionFuncs as dcf

# -------------------------------------------------------------#
# OPTIONS
# -------------------------------------------------------------#

# Set path to folder containing full EddyPro cospectra
path_cospec = '/Path/To/Cospectra/';

# Set path to file containing full half-hour output (contains RH)
path_fullout = '/Path/To/Fullout/file.csv';

# Name of sensor (used for saving file)
sensorName = 'SaveSensorName'

# Select Gases as list of strings.
# Current options 'co2','h2o'
# No need to include 'ts' because it is added later
gases = ['co2','h2o']

# Select relative humidity bin type
# Current options are 'none' and 'quantile', but initialization below
# could easily be altered to construct different bins
# List of strings, same size as gases.
bin_opt = ['none','quantile']

# Set start and end dates ('YYYY-MM-DD')
start_date = '2015-07-05'  # Inclusive: starts at 00:00 of given day
end_date = '2015-10-05'    # Inclusive: ends at 00:00 of the following day

# File to save out cospectra, because loading takes several minutes
# If new date range is selected, cospectra will be re-loaded and saved separately

npzout = '{}_DirectCospec_{}_{}'.format(sensorName,start_date,end_date)

# -------------------------------------------------------------#
# SANITY CHECKS
# -------------------------------------------------------------#
if len(gases)!=len(bin_opt):
	print('Must choose one bin option per selected gas!')

for i in range(len(gases)):
	if gases[i] not in ['co2','h2o','ch4']:
		print('Selected gases must be co2, h2o, or ch4!')
	if bin_opt[i] not in ['none','quantile']:
		print('Selected bin_opt must be iether none or quantile!')

if datetime.strptime(end_date,'%Y-%m-%d') < datetime.strptime(start_date,'%Y-%m-%d'):
	print('End date must be after start date!')
	
# -------------------------------------------------------------#
# RUN SPECTRAL CALCULATIONS
# -------------------------------------------------------------#

# Variables to load: gases and ts
vars = gases + ['ts']

# Load cospectra. This can be slow depending on the time range.
# Often requires quite a bit of RAM since all cospectra are loaded simultaneously
# Returns co2, h2o, ts cospectra and the corresponding frequencies

allfiles = os.listdir('./')
if npzout+'.npz' in allfiles:
	loaded=np.load(npzout+'.npz')
	cospec=loaded['cospec']
	time=loaded['time']
	freq=loaded['freq']
	print('Loaded cospectra from {}'.format(npzout+'.npz'))
else:
	cospec,time,freq = dcf.Load_EP_Fullcospectra(path_cospec,start_date,end_date,vars)
	np.savez(npzout,cospec=cospec,time=time,freq=freq)
	print('Saved cospectra to {}'.format(npzout+'.npz'))

bindex,btime,data,header = dcf.make_RH_bins(path_fullout,start_date,end_date,bin_opt)

# Insert index here to apply filters
# Set points to NaN that should not be used (e.g. cospec[:,~index,:] = np.nan)
# Calculate Direct spectral correction

cospec_indexed = cospec*1.

# Example of wind direction index:
#WDi = np.array([True if hstr=='wind_dir' else False for hstr in header])
#index = np.abs(data[:,WDi]-180)<60
#index = np.ndarray.flatten(index)
#cospec_indexed[:,~index,:]=np.nan

direct_scf = dcf.Calculate_Direct_SCF(cospec_indexed,freq,bindex)

# -------------------------------------------------------------#
# PLOT
# -------------------------------------------------------------#
plt.figure(1)
plt.clf
for i in range(len(vars)):
	plt.semilogx(freq, np.nanmedian(cospec_indexed[:,:,i],1),label=vars[i])
plt.legend(vars,fontsize=14)

plt.figure(2)
plt.clf
plt.plot(time,direct_scf,'.')
plt.ylabel('Direct Correction',fontsize=16)
plt.legend(gases)

plt.show()

	
