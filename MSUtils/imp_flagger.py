import matplotlib.pyplot as plt
import numpy as np
import itertools 
import time as T
import pandas as pd
import astropy
from astropy.time import Time
#import katdal as kd
from pyrap.tables import table
import collections
import matplotlib


def plot_gain_table(gain_table, plt_scale, plot_file):
#Reads in the gain table; makes plots of gain amplitude/phase with 
#time/frequency for all the antennas.

#Read in the table
   G_tab = table(gain_table)
   print 'Reading gain table:', gain_table, '\n'

#Get number of antennas (needed for plotting)       
   Ant_n1 = G_tab.getcol("ANTENNA1")
   Ant_n2 = G_tab.getcol("ANTENNA2")
   Ant_list = list(set(np.append(Ant_n1,Ant_n2)))
   N_ants = len(Ant_list)
   print 'Number of Antennas to plot:', N_ants    
   
#Read in the flags
   flags = G_tab.getcol("FLAG")
   
#Read in the solutions
   Gsols = G_tab.getcol("CPARAM")
   nchans = Gsols.shape[1]
   ntint  = Gsols.shape[0]
   npol = Gsols.shape[2]        
#Read in the error.
   Gsols_err = G_tab.getcol("PARAMERR")

#Read in the timestamps
 
   Gsols_time = G_tab.getcol("TIME")

#Prepare for plotting; store gain amp and phases.
 
   Gsols_HH = Gsols[:,:,0]
   Gsols_VV = Gsols[:,:,1]
   Gsols_HH_flag = np.ma.masked_array(Gsols_HH, mask=flags[:,:,0])
   Gsols_VV_flag = np.ma.masked_array(Gsols_VV, mask=flags[:,:,1])
   Gsols_HH_amp = abs(Gsols_HH_flag)
   Gsols_VV_amp = abs(Gsols_VV_flag)
   Gsols_HH_ph = np.angle(Gsols_HH_flag,deg=True)
   Gsols_VV_ph = np.angle(Gsols_VV_flag,deg=True)
   


#Plotting
#Plot in a more or less square grid.
   nplts = int(np.sqrt(N_ants))       #(if non zero remainder, add one)

#Set Global matplotlib options
   matplotlib.rcParams['lines.markersize'] = 10.0
   matplotlib.rcParams['xtick.major.size'] = 11.5
   matplotlib.rcParams['ytick.major.size'] = 11.5
   matplotlib.rcParams['ytick.direction'] = 'out'
   matplotlib.rcParams['xtick.direction'] = 'out'
   matplotlib.rcParams['font.size'] = 11.0

   f, axarr = plt.subplots(nplts, nplts, dpi=800, figsize=(nplts*plt_scale,nplts*plt_scale))
#Plot amplitudes first
   for ant in xrange(N_ants):
       axarr[ant // nplts, ant % nplts].plot(Gsols_HH_amp[ant],'g.')
       axarr[ant // nplts, ant % nplts].plot(Gsols_VV_amp[ant],'r.')
       axarr[ant // nplts, ant % nplts].set_title('Antenna '+str(ant))
       axarr[ant // nplts, ant % nplts].set_ylim(np.round(np.min(Gsols_HH_amp),1)-0.05,np.round(np.max(Gsols_HH_amp),1)+0.05)
#   f.savefig(output_dir + "/%s-AUTOCORR-FIELD-%s-CORR-%d.png" %
#            (os.path.basename(ms),
#            source_names[field_id],
#            corr))
   f.savefig(plot_file) 
   plt.close(f)

   return

plot_gain_table('meerkathi-12A-405.sb11468652.eb11592823.56163.7871827662-1gc1.B0', 6, 'temp_plot.png')
