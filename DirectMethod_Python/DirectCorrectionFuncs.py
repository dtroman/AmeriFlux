# DirectCorrectionFuncs.py
# Contains functions used by DirectCorrectionMaster.py

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import os
import csv
import re


def Calculate_Direct_SCF(cospectra, freq, bindex):
    Nfreq = cospectra.shape[0]
    Nt = cospectra.shape[1]
    Nvars = bindex.shape[1]
    
    # Number of frequency bins. 45 because it is the number of binned 
    # frequencies in the EddyPro binned cospectra (at least for 30 min 
    # averaging)
    Nfb = 45
    scf = np.zeros_like(bindex)*np.nan
    
    # Figure out points to bin around (similar to 'binned frequency')
    # frequencies in the EddyPro binned cospectra (at least for 30 min
    # averaging)
    tmpi = np.logspace(np.log10(1),np.log10(Nfreq),Nfb) 
    i=np.array(np.round(tmpi),dtype=int)
    def idx(v):
        return np.array(np.hstack((1,np.diff(v))),dtype=int)
    while any(idx(i)==0):
        fi = np.where(idx(i)==0)[0]
        i[fi[0]:np.int(fi[0]+np.floor(Nfb/2))+1]=i[fi[0]:int(fi[0]+np.floor(Nfb/2))+1]+1
    
    for vi in range(Nvars):
        bindex_cur = bindex[:,vi]
        # (Nbins below) Minus one because 0 is for points that do not fall in a bin
        Nbins = len(np.unique(bindex_cur))-1; 
        for b in range(1,Nbins+1):
            bi = bindex_cur==b
        
            gasspec = np.nanmedian(cospectra[:,bi,vi],1)
            Hspec   = np.nanmedian(cospectra[:,bi,Nvars],1)
            
            fbinmin = np.hstack((0,i[:-1]+np.diff(i)/2))
            fbinmax = np.hstack((i[:-1]+np.diff(i)/2,Nfreq+1))
            allidx = np.arange(1,Nfreq+1)
            
            gasspec_bin = np.zeros(Nfb) * np.nan
            Hspec_bin = np.zeros(Nfb) * np.nan
            freq_bin = np.zeros(Nfb) * np.nan
            
            for fb in range(len(fbinmin)):
                fbi = (allidx>=fbinmin[fb]) & (allidx<fbinmax[fb])
                freq_bin[fb] = np.median(freq[fbi])
                gasspec_bin[fb] = np.mean(gasspec[fbi])
                Hspec_bin[fb] = np.mean(Hspec[fbi])
                
            # Calculate transfer functions
            tfbin = gasspec_bin/Hspec_bin  # Fratini et al., 2012 Eq. 1
            
            # Interpolate transfer function back to high frequencies
            tfinterp = np.interp(np.arange(1,Nfreq+1),i,tfbin)
            # Limit analysis to high frequencies and remove negative points
            # If there is aliasing in the cospectra, this is where changes are required
            lowtf = tfinterp<0
            if np.sum(lowtf>0)>0:
                tfinterp[tfinterp<0] = np.min(tfinterp[tfinterp>0])
            if np.sum(tfinterp>=1)>0:
                lowest = np.where(tfinterp>=1)[0][-1]
                tfinterp[:(lowest+1)]=1
            
            tfinterp_stack = np.swapaxes(np.tile(tfinterp,(Nt,1)),0,1)
            freq_stack = np.swapaxes(np.tile(freq,(Nt,1)),0,1)
            
            origspec = cospectra[:,:,Nvars]
            
            df = freq[1]-freq[0]
            # Equation 3 of Fratini et al., 2012 (no square root)
            Fl_numerator   = np.sum(origspec*df/freq_stack,0)
            Fl_denominator = np.sum(origspec*df/freq_stack*tfinterp_stack,0)
            
            scf_all = Fl_numerator/Fl_denominator
            scf[bi,vi] = scf_all[bi]

    scf[scf<0]  = np.nan
    scf[scf>50] = np.nan
    
    return scf
    
def make_RH_bins(path,start_date,end_date,bin_opt):
    """
    make_RH_bins(path,start_date,end_date,bin_opt) takes full EddyPro outputs
    and returns a vector of bin numbers.
    path is the path and file name of the EddyPro fulloutput file
    start_date and end_date have the format 'YYYY-MM-DD'
    bin_opt is 'none' or 'quantile' - could be adapted to include more
    """
    
    # Load full output to read relative humidities (RH)
    time=[]
    with open(path, "r") as f:
        reader = csv.reader(f,delimiter=',')
        ct=1
        for row in reader:
            if ct==2:
                header = row
            elif ct>3:
                curtime = datetime.strptime('{} {}'.format(row[1],row[2]),
                                            '%Y-%m-%d %H:%M')
                time.append(curtime)
            ct+=1
    
    # Remove text columns from data and corresponding headers
    header = header[3:]
    data = np.genfromtxt(path,delimiter=',',skip_header=3)  
    data = data[:,3:]
    
    sdatetime = datetime.strptime(start_date,'%Y-%m-%d')
    edatetime = datetime.strptime(end_date,'%Y-%m-%d')
    edatetime = edatetime+timedelta(1) # Add a day to include whole end_date
    # Cut to selected dates
    timei = [(t>sdatetime)&(t<=edatetime) for t in time]
    data = data[timei,:]
    
    RHi = [h=='RH' for h in header]
    RHi = np.where(RHi)[0][0]
    RH = np.ndarray.flatten(data[:,RHi])
    nni = ~np.isnan(RH)
    # Get the indices for each selected option
    bindex = np.zeros((len(RH),len(bin_opt)))
    for boi in range(len(bin_opt)):
        # Apply RH bin options to get bin limits
        if bin_opt[boi].lower()=='quantile':
            RH_bins = np.hstack((np.nanmin(RH)-0.01,
                                np.quantile(RH[nni],[0.25,0.5,0.75]),
                                np.nanmax(RH)+0.01))
        elif bin_opt[boi].lower()=='none':
            RH_bins = np.array([0,100])
    
        # Label each point with a bin (e.g. 0=exclude, 1=first bin, etc.)
        for bi in range(1,len(RH_bins)):
            thisbin = (RH>=RH_bins[bi-1]) & (RH<RH_bins[bi])
            bindex[thisbin,boi] = bi
    
    return bindex, time, data, header 


def Load_EP_Fullcospectra(path,start_day,end_day,variable):
    """
    LoadDavisSpectra loads the Ameriflux Davis spectral outputs from EddyPro
    Start Day and End Day are input as date strings ('yyyy-mm-dd')
    Path is the path from the current folder to get to folder with all the
    data (organized by date)
        Options: 'spec','cospec','both'    
    """
    
    # Number of days selected
    sday = datetime.strptime(start_day,'%Y-%m-%d')
    eday = datetime.strptime(end_day,'%Y-%m-%d')
    Nday = (eday-sday).days +1
    
    if Nday <= 0:
        print('WARNING!! End day is before start day!')
        
    Nvars = len(variable)

    allf = os.listdir(path)
    fnames = [f for f in allf if f.endswith('.csv')]
    
    # Read first file to get info (meta)    
    spec, timeseries, header, meta1 = read_cospectrum(path,[fnames[0]])
    Hz = meta1[0]
    avg_period = meta1[3]
    nseg = np.int(24*60/avg_period)
    ppf = np.int(2**np.floor(np.log2(avg_period*60*Hz/2)))

    df = Hz/2/ppf
    freq = np.arange(df,Hz/2+df,df)
    
    # spec shape: [frequency,time,variables]
    spec=np.zeros((ppf,np.int(Nday*(24*60/avg_period)),Nvars))*np.nan
    spec_time=[]

    tct = -1 # Time counter
    for d in range(Nday):
        for h in range(nseg):
                        tct+=1
                        curtime = sday+timedelta(d,0,0,0,avg_period*(h+1))
                        spec_time.append(curtime)
                        hstr = (curtime).strftime('%H%M')

                        daystr = curtime.strftime('%Y-%m-%d')
                        daystr2 = curtime.strftime('%Y%m%d')
                        print('Loading... {} {}'.format(daystr,hstr))

                        # See if file exists
                        matchi = np.array(['{}-{}'.format(daystr2,hstr) in f for f in fnames])

                        if np.sum(matchi)>0:
                                matchi = np.where(matchi)[0][0]
                                spec_day, spec_time_day, header_day, meta_day = read_cospectrum(path,[fnames[matchi]])
                                spec_day = spec_day[0]

                                for vi in range(Nvars):
                                        gasheader = 'f_nat*cospec(w_{})'.format(variable[vi])
                                        vmatchi = np.array([gasheader in h for h in header_day])
                                        if np.sum(vmatchi)>0:
                                                vmatchi = np.where(vmatchi)[0][0]
                                                spec[:,tct,vi] = spec_day[:,vmatchi]

                                        else:
                                                print('And there was a problem!')                   
    
    return spec, spec_time, freq
    
def read_cospectrum(path,d):
    """
    batch loader for EddyPro binned (co)spectra results
    inputs: 
    path = directory where files are located
    d = list of files names to load
    Assumes that there is a date in the file name with the format 
        yyyymmdd-HHMM
    
    outputs:
    spec = stacked output of all frequency bins x columns of (co)spectra
    timeseries = list of datetime objects of each file used
    header = list of column names
    meta = [sampling frequency,sample height, wind speed, averaging period]
    """
    spec = []
    timeseries = []
    for i in range(len(d)):
        filename = path + d[i]

        with open(filename, "r") as f:
            reader = csv.reader(f,delimiter=',')
            ct=1
            for row in reader:
                if ct==6:
                    Hz = float(row[0].split('_')[-1])
                elif ct==7:
                    height = float(row[0].split('_')[-1])
                elif ct==8:
                    ws = float(row[0].split('_')[-1])
                elif ct==9:
                    avg_period = float(row[0].split('_')[-1])
                elif ct==13:
                    header = row
                elif ct>13:
                    break
                ct+=1
                
        meta = [Hz,height,ws,avg_period]
            
        thisspec = np.genfromtxt(filename,delimiter=',',skip_header=13)
        spec.append(thisspec)
        thistime = re.findall('\d{8}-\d{4}',filename)[0]
        thisdate = datetime.strptime(thistime,'%Y%m%d-%H%M')
        timeseries.append(thisdate)                         
                
    return spec, timeseries, header, meta
    
    
