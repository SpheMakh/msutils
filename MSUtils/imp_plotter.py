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


def plot_bandpass_table(gain_table, plt_scale, plt_dpi, plot_file):
#Reads in the gain table; makes plots of gain amplitude/phase with 
#time/frequency for all the antennas.

#Read in the table
   G_tab = table(gain_table)
   print 'Reading gain table:', gain_table, '\n'
   G_tab_names = table(gain_table+"::ANTENNA")
   ant_names = G_tab_names.getcol("NAME")
   print "Antennas present in the table:", ant_names

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
   nscans  = Gsols.shape[0]/N_ants
   npol = Gsols.shape[2]        
#Read in the error.
   Gsols_err = G_tab.getcol("PARAMERR")

#Read in the timestamps
 
   Gsols_time = G_tab.getcol("TIME")

#Prepare for plotting; store gain amp and phases.
 
   Gsols_HH = Gsols[:,:,0]
   Gsols_VV = Gsols[:,:,1]
   print "Shape of solutions:", np.shape(Gsols_HH), np.shape(Gsols_VV)
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
   matplotlib.rcParams['lines.markersize'] = 4.0
   matplotlib.rcParams['xtick.major.size'] = 5.0
   matplotlib.rcParams['ytick.major.size'] = 5.0
   matplotlib.rcParams['ytick.direction'] = 'out'
   matplotlib.rcParams['xtick.direction'] = 'out'
   matplotlib.rcParams['font.size'] = 15.0
   matplotlib.rcParams['axes.linewidth'] = 3.0
   matplotlib.rcParams['legend.framealpha'] = 0.5
#   matplotlib.rcParams[]
   
   plt_dpi=800
   f, axarr = plt.subplots(nplts, nplts, dpi=plt_dpi, figsize=(nplts*plt_scale,nplts*plt_scale))
   f.text(0.5,0.94,"Bandpass Plot",ha='center',fontsize=40)
   f.text(0.5, 0.04, 'Channel', ha='center',fontsize=30)
   f.text(0.04, 0.5, 'Gain Amp', va='center', rotation='vertical',fontsize=30)

#Plot amplitudes first
   for ant in xrange(N_ants):
       for scan in xrange(nscans):
           axarr[ant // nplts, ant % nplts].plot(Gsols_HH_amp[ant+N_ants*scan],'g.', label='CORR0')
           axarr[ant // nplts, ant % nplts].plot(Gsols_VV_amp[ant+N_ants*scan],'r.', label='CORR1')
           if scan==0:
              axarr[ant // nplts, ant % nplts].legend()
       axarr[ant // nplts, ant % nplts].set_title('Antenna '+str(ant_names[ant]))
       axarr[ant // nplts, ant % nplts].set_ylim(np.round(np.min(Gsols_HH_amp),1)-0.05,np.round(np.max(Gsols_HH_amp),1)+0.05)
   plot_f = plot_file+"bandpass-amp.png"
   f.savefig(plot_f) 
   plt.close(f)

   f, axarr = plt.subplots(nplts, nplts, dpi=plt_dpi, figsize=(nplts*plt_scale,nplts*plt_scale))
   f.text(0.5,0.94,"Bandpass Plot",ha='center',fontsize=40)
   f.text(0.5, 0.04, 'Channel', ha='center',fontsize=30)
   f.text(0.04, 0.5, 'Gain Phase', va='center', rotation='vertical',fontsize=30)

#Plot phases now
   for ant in xrange(N_ants):
       for scan in xrange(nscans):
           axarr[ant // nplts, ant % nplts].plot(Gsols_HH_ph[ant+N_ants*scan],'g.', label='CORR0')
           axarr[ant // nplts, ant % nplts].plot(Gsols_VV_ph[ant+N_ants*scan],'r.', label='CORR1')
           if scan==0:
              axarr[ant // nplts, ant % nplts].legend()
       axarr[ant // nplts, ant % nplts].set_title('Antenna '+str(ant))
       axarr[ant // nplts, ant % nplts].set_ylim(np.round(np.min(Gsols_HH_ph),5)-0.5,np.round(np.max(Gsols_HH_ph),5)+0.5)
   plot_f = plot_file+"-bandpass-phase.png"
   f.savefig(plot_f)
   plt.close(f)

   return

def plot_gain_table(gain_table, plt_scale, plt_dpi, plot_file):
#Reads in the gain table; makes plots of gain amplitude/phase with 
#time/frequency for all the antennas.

#Read in the table
   G_tab = table(gain_table)
   print 'Reading gain table:', gain_table, '\n'
   G_tab_names = table(gain_table+"::ANTENNA")
   ant_names = G_tab_names.getcol("NAME")
   print "Antennas present in the table:", ant_names

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
   nscans  = Gsols.shape[0]/N_ants
   print "Number of scans:", nscans
   npol = Gsols.shape[2]
#Read in the error.
   Gsols_err = G_tab.getcol("PARAMERR")

#Read in the timestamps

   Gsols_time = G_tab.getcol("TIME")

#Prepare for plotting; store gain amp and phases.

   Gsols_HH = Gsols[:,:,0]
   Gsols_VV = Gsols[:,:,1]
   print "Shape of solutions:", np.shape(Gsols_HH), np.shape(Gsols_VV)
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
   matplotlib.rcParams['lines.markersize'] = 15.0
   matplotlib.rcParams['xtick.major.size'] = 5.0
   matplotlib.rcParams['ytick.major.size'] = 5.0
   matplotlib.rcParams['ytick.direction'] = 'out'
   matplotlib.rcParams['xtick.direction'] = 'out'
   matplotlib.rcParams['font.size'] = 15.0
   matplotlib.rcParams['axes.linewidth'] = 3.0
   matplotlib.rcParams['legend.framealpha'] = 0.5
  
#Make a single plot for gain amplitude and phase, colourize different fields?


   obs_time = G_tab.getcol("TIME")

#Scale to start time and scale to hours
   obs_time = obs_time-obs_time[0]

   obs_time = obs_time/3600.0







   f, axarr = plt.subplots(nplts, nplts, dpi=plt_dpi, figsize=(nplts*plt_scale,nplts*plt_scale))
#Plot amplitudes first
   f.text(0.5,0.94,"Gain Amplitude Plot",ha='center',fontsize=40)
   f.text(0.5, 0.04, 'Time(Hours)', ha='center',fontsize=30)
   f.text(0.04, 0.5, 'Gain Amp', va='center', rotation='vertical',fontsize=30)
   for ant in xrange(N_ants):
       for scan in xrange(nscans):
           axarr[ant // nplts, ant % nplts].plot(obs_time[ant+N_ants*scan], Gsols_HH_amp[ant+N_ants*scan],'g.', label='CORR0')
           axarr[ant // nplts, ant % nplts].plot(obs_time[ant+N_ants*scan], Gsols_VV_amp[ant+N_ants*scan],'r.', label='CORR1')
           if scan==0:
               axarr[ant // nplts, ant % nplts].legend()
           axarr[ant // nplts, ant % nplts].set_title('Antenna '+str(ant_names[ant]))
      #     axarr[ant // nplts, ant % nplts].set_xlabel('Time (Hours)')
      #     axarr[ant // nplts, ant % nplts].set_ylabel('Gain Amplitude')
           axarr[ant // nplts, ant % nplts].set_ylim(0,20)
#       axarr[ant // nplts, ant % nplts].set_ylim(np.round(np.min(Gsols_HH_amp),1)-0.05,np.round(np.max(Gsols_HH_amp),1)+0.05)
   plot_f = plot_file+"-gain-amp.png"
   f.savefig(plot_f)
   plt.close(f)

   f, axarr = plt.subplots(nplts, nplts, dpi=plt_dpi, figsize=(nplts*plt_scale,nplts*plt_scale))
   #Plot phases first
   f.text(0.5,0.94,"Gain Phase Plot",ha='center',fontsize=40)
   f.text(0.5, 0.04, 'Time(Hours)', ha='center',fontsize=30)
   f.text(0.04, 0.5, 'Gain Phase(Degrees)', va='center', rotation='vertical',fontsize=30)
   for ant in xrange(N_ants):
       for scan in xrange(nscans):
           axarr[ant // nplts, ant % nplts].plot(obs_time[ant+N_ants*scan], Gsols_HH_ph[ant+N_ants*scan],'g.', label='CORR0')
           axarr[ant // nplts, ant % nplts].plot(obs_time[ant+N_ants*scan], Gsols_VV_ph[ant+N_ants*scan],'r.', label='CORR1')
           if scan==0:
               axarr[ant // nplts, ant % nplts].legend()
           axarr[ant // nplts, ant % nplts].set_title('Antenna '+str(ant_names[ant]))
      #     axarr[ant // nplts, ant % nplts].set_xlabel('Time (Hours)')
      #     axarr[ant // nplts, ant % nplts].set_ylabel('Gain Amplitude')
      #     axarr[ant // nplts, ant % nplts].set_ylim(0,20)
           axarr[ant // nplts, ant % nplts].set_ylim(np.round(np.min(Gsols_HH_ph),5)-0.5,np.round(np.max(Gsols_HH_ph),5)+0.5)
   plot_f = plot_file+"-gain-ph.png"
   f.savefig(plot_f)
   plt.close(f)

   return

def plot_delay_table(gain_table, plt_scale, plt_dpi, plot_file):
#Reads in the gain table; makes plots of gain amplitude/phase with 
#time/frequency for all the antennas.

#Read in the table
   G_tab = table(gain_table)
   print 'Reading gain table:', gain_table, '\n'
   G_tab_names = table(gain_table+"::ANTENNA")
   ant_names = G_tab_names.getcol("NAME")
   print "Antennas present in the table:", ant_names

#Get number of antennas (needed for plotting)       
   Ant_n1 = G_tab.getcol("ANTENNA1")
   Ant_n2 = G_tab.getcol("ANTENNA2")
  # Ant_list = list(set(np.append(Ant_n1,Ant_n2)))
   N_ants = len(ant_names)
   print 'Number of Antennas to plot:', N_ants

#Read in the flags
   flags = G_tab.getcol("FLAG")

#Read in the solutions
   Gsols = G_tab.getcol("FPARAM")
   nchans = Gsols.shape[1]
   nscans  = Gsols.shape[0]/(N_ants+1)   #+1 since one of the antennas is the reference.
   print "Number of scans:", nscans
   npol = Gsols.shape[2]
#Read in the error.
   Gsols_err = G_tab.getcol("PARAMERR")

#Read in the timestamps

   Gsols_time = G_tab.getcol("TIME")

#Prepare for plotting; store gain amp and phases.

   Gsols_HH = Gsols[:,:,0]
   Gsols_VV = Gsols[:,:,1]
   print "Shape of solutions:", np.shape(Gsols_HH), np.shape(Gsols_VV)
   Gsols_HH_plt = np.ma.masked_array(Gsols_HH, mask=flags[:,:,0])
   Gsols_VV_plt = np.ma.masked_array(Gsols_VV, mask=flags[:,:,1])
   #Plotting
#Plot in a more or less square grid.
   nplts = int(np.sqrt(N_ants+1))       #(if non zero remainder, add one)
   print "Number of plots:", nplts*nplts

#Set Global matplotlib options
   matplotlib.rcParams['lines.markersize'] = 15.0
   matplotlib.rcParams['xtick.major.size'] = 5.0
   matplotlib.rcParams['ytick.major.size'] = 5.0
   matplotlib.rcParams['ytick.direction'] = 'out'
   matplotlib.rcParams['xtick.direction'] = 'out'
   matplotlib.rcParams['font.size'] = 15.0
   matplotlib.rcParams['axes.linewidth'] = 3.0
   matplotlib.rcParams['legend.framealpha'] = 0.5
#   matplotlib.rcParams['legend.frameon']='false'
#Make a single plot for gain amplitude and phase, colourize different fields?


   obs_time = G_tab.getcol("TIME")

#Scale to start time and scale to hours
   obs_time = obs_time-obs_time[0]

   obs_time = obs_time/3600.0
   
   f, axarr = plt.subplots(nplts, nplts, dpi=plt_dpi, figsize=(nplts*plt_scale,nplts*plt_scale)) 
   
   f.text(0.5,0.94,"Delay Calibration Plot",ha='center',fontsize=40)
   f.text(0.5, 0.04, 'Time (Hours)', ha='center',fontsize=30)
   f.text(0.04, 0.5, 'Delay (ns)', va='center', rotation='vertical',fontsize=30)
   for ant in xrange(N_ants):
       for scan in xrange(nscans):
           axarr[ant // nplts, ant % nplts].plot(obs_time[ant+N_ants*scan], Gsols_HH_plt[ant+N_ants*scan],'g.', label='CORR0')
           axarr[ant // nplts, ant % nplts].plot(obs_time[ant+N_ants*scan], Gsols_VV_plt[ant+N_ants*scan],'r.', label='CORR1')
           if scan==0:
               axarr[ant // nplts, ant % nplts].legend()
           axarr[ant // nplts, ant % nplts].set_title('Antenna '+str(ant_names[ant]))
      #     axarr[ant // nplts, ant % nplts].set_xlabel('Time (Hours)')
      #     axarr[ant // nplts, ant % nplts].set_ylabel('Gain Amplitude')
      #     axarr[ant // nplts, ant % nplts].set_ylim(0,20)
           axarr[ant // nplts, ant % nplts].set_ylim(np.round(np.min(Gsols_HH_plt),5)-0.5,np.round(np.max(Gsols_HH_plt),5)+0.5)
   plot_f = plot_file+"-delay.png"
   f.savefig(plot_f)
   plt.close(f)

   return

def gain_plotter(caltable,typ,outfile,plt_scale=6,plt_dpi=600):
#"Plots gaintables. Types should be 'delay', 'gain' or 'bandpass'. Scale is the plot scale, dpi is the plot dpi"
    if (typ=='delay'):
        plot_delay_table(caltable,plt_scale,plt_dpi,outfile)
    if (typ=='gain'):
        plot_gain_table(caltable,plt_scale,plt_dpi,outfile)
    if (typ=='bandpass'):
        plot_bandpass_table(caltable,plt_scale,plt_dpi,outfile)

    return

gain_plotter('meerkathi-12A-405.sb11468652.eb11592823.56163.7871827662-1gc1.K0','delay','fin_test')
#plot_bandpass_table('meerkathi-12A-405.sb11468652.eb11592823.56163.7871827662-1gc1.B0', 6, 800,'temp_plot')
#plot_gain_table('meerkathi-12A-405.sb11468652.eb11592823.56163.7871827662-1gc1.F0', 6, 'temp_plot')
#plot_delay_table('meerkathi-12A-405.sb11468652.eb11592823.56163.7871827662-1gc1.K0', 6, 'temp_plot')
