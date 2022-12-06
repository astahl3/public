0#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 16:28:13 2022

@author: astahl3
"""
import re
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import os
import pathlib
import math

####################################################
######## ENTER DESIRED NUMERICAL VARIABLES #########
####################################################

# 0 = display name
# 1 = pattern match for log file (copy-paste entire process_0.log line up until the desired number)
var_01 =  ['loop time (sec)',               r' - Loop time  \(mpi   time\):.*?([0-9.-]+)']
var_02 =  ['fastest particle (level 0)',    r'   ->Fastest particle \(int\) in Level 0:.*?([0-9.-]+)']
var_03 =  ['fastest particle (level 1)',    r'   ->Fastest particle \(int\) in Level 1:.*?([0-9.-]+)']
var_04 =  ['fastest particle (level 2)',    r'   ->Fastest particle \(int\) in Level 2:.*?([0-9.-]+)']
var_05 =  ['fastest particle (level 3)',    r'   ->Fastest particle \(int\) in Level 3:.*?([0-9.-]+)']

variables  =    [
                var_01,
                #var_02,
                #var_03,
                #var_04,
                #var_05,
                ]
####################################################
############## ENTER PATH TO LOG FILE ##############
####################################################

#log_path = '/home/astahl3/aikef_output/Juno/Juno_n10_r0_dt10/bin/process_0.log'
#log_path = '/home/astahl3/aikef_output/Juno/Juno_n10_r1_dt10_Duling/bin/process_0.log'
#log_path = '/home/astahl3/aikef_output/Juno/Juno_n14_r3_dtplus/bin/process_0.log'
log_path = '/home/astahl3/aikef_output/Juno/Juno_n12_r1_dt10/bin/process_0.log'
#log_path = '/data1/astahl/aikef_output/G2/G2_master_r0/bin/process_0.log'
#log_path = '/home/astahl3/aikef_output/Juno/Juno_n10_r1_dt10_D2/bin/process_0.log'

####################################################
########## EXTRACT BASIC INFO ABOUT RUNS ###########
####################################################
TL_per_run = []
dt_per_run = []
TL_save = []
TL_save.append(0)
run_count = 1
TL_temp = 0
total_lines = 0
cnt = 0

###### Saving basic info about number of runs, timesteps of each run in log file, starting save position
with open(log_path) as log_file:
    pattern_match_dt = re.compile(r' dt = .*?([0-9.-]+)')            
    pattern_match_time = re.compile(r' STATISTICS OF TL .*?([0-9.-]+)')
    pattern_match_save = re.compile(r' Continuing at TL .*?([0-9.-]+)')

    for line in log_file:
        total_lines = total_lines + 1 # keep count of total number of lines read
        cnt = cnt + 1 # dummy counter for print statements
        match_dt = pattern_match_dt.match(line) # extract the AIKEF timestep "dt" used in current run
        match_time = pattern_match_time.match(line) # extract timestep n
        match_save = pattern_match_save.match(line) # extract timestep for restart from state
        
        if match_time:
            TL_temp = match_time.group(1) # save most recent timestep in current run
        
        # If true, this means new run is starting
        # Save the final timestep "TL" of previous run
        # Save the AIKEF timestep "dt" used in current run
        if match_dt:
            dt_per_run.append(match_dt.group(1)) # extract dt used for current run
            TL_per_run.append(TL_temp) # save the associated (i.e., max timestep of prior run)
        if match_save:
            TL_save.append(match_save.group(1))
            run_count = run_count + 1
            
        if cnt == 100000:
            print(total_lines, 'lines scanned during run counter')
            cnt = 0
            
    TL_per_run.append(TL_temp) # store total number of timesteps achieved in each run
    TL_per_run.pop(0) # remove first entry as it doesn't corresond to any run

TL_save = np.float_(TL_save)
TL_save = np.int_(TL_save)
TL_per_run = np.float_(TL_per_run)
TL_per_run = np.int_(TL_per_run)

###### Print summary of each run included in the log file
print('\n')
for k in range(0,run_count): 
    TL_per_run[k] = TL_per_run[k] - TL_save[k] 
    print('run #', k+1, 'used dt =', dt_per_run[k], 'totalling ', TL_per_run[k], 'timesteps')

print('\n')
print('TOTAL LINES IN LOG FILE: ', total_lines, '\n')

####################################################
########### STRINGS TO MATCH IN LOG FILE ###########
####################################################

    ##############################################################################################      
    # Extracts the number (as a float) appearing after a particular line                         #
    # The string of the line must be copied exactly, including leading spaces, from              #
    # start of line until some point before the desired numerical value                          #
    # Use backslash before parenthesis, e.g.:                                                    #   
    #    if desired line is:  "   -> vel (int) in level 0 = " (wihthout the quotes)              #
    #    enter this in  re.compile():  "    -> vel \(int\) in level 0 = " (without the quotes)   #
    ##############################################################################################

pattern_start = re.compile(variables[0][1])
                   
###### Save the desired variable or value at each timestep
save_var = []
with open(log_path) as log_file:
    k = 0
    cnt = 0
    for line in log_file:
        match_var = pattern_start.match(line)
        if match_var:
            #print('match.group(1)')
            save_var.append(match_var.group(1))
            #print(match.group(1))
            #result.append(match.group(1))
            k = k + 1
            cnt = cnt + 1
        if k == 500:
            print(cnt, ' timesteps completed during variable extraction')
            k = 0
               
###### PLOT ALL TIMESTEPS ######
var_plot = np.float_(save_var)            
fig, ax = plt.subplots(1,1,sharex=True)
plt.subplots_adjust(hspace = 0.05)
ax.plot(np.linspace(1,cnt,cnt,endpoint=True),var_plot)
ax.set_yticks(np.round(np.linspace(1,np.max(var_plot),10,endpoint=True)))
ax.set_ylabel('maximum normalized v')
ax.set_xlabel('timestep')
ax.set_title('all timesteps')
#plt.savefig('/home/astahl3/plots/Juno_n10_r1_dt10_Duling_vel.png')

###### MAKE PLOTS FOR EACH RUN ######

for k in range(0,run_count):
    fig, ax = plt.subplots(1,1,sharex=True)
    plt.subplots_adjust(hspace = 0.05)
    temp_plot = var_plot[TL_save[k]:TL_save[k]+TL_per_run[k]:1]
    ax.plot(np.linspace(TL_save[k],TL_per_run[k],TL_per_run[k],endpoint=True),temp_plot)
    ax.set_title('run # '+str(k+1))
    ax.set_xlabel('timestep')
    ax.set_ylabel(variables[0][0])   
    
    
###### PLOT TIMESTEPS LAST RUN ONLY ######
#TL_per_run = np.float_(TL_per_run)
#num_runs = len(TL_per_run)
#TL_last_end = np.sum(TL_per_run)
#TL_last_start = np.sum(TL_per_run[0:num_runs-1])
#steps = round(TL_last_end)-round(TL_last_start)
#var_last_plot = var_plot[round(TL_last_start):steps:1]

#fig, ax = plt.subplots(1,1,sharex=True)
#plt.subplots_adjust(hspace = 0.05)
#ax.plot(np.linspace(1,steps,steps,endpoint=True),var_last_plot)
#ax.set_yticks(np.round(np.linspace(1,np.max(var_plot),10,endpoint=True)))
#ax.set_ylabel(variables[0][0])
#ax.set_xlabel('timestep')

