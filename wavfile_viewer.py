#-------------------------------------------------------------------------------
# Name:        wavfile viewer
# Purpose:
#
# Author:      kdhollan
#
# Created:     28-05-2013
# Copyright:   (c) kdhollan 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# GIT password = FEICompany2013


import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import os
import subprocess

from scipy.io import wavfile
from scipy import stats
from pylab import *














# Definition of global parameters
tracelist =    ['Copy of input signal',
                'HPD output (including datapath delay)',
                'Bilateral filter output',
                'Bilateral filter error signal (input for threshold detector)',
                'FIR1 output (1st order)',
                'FIR2 output (3rd order)',
                'Bilat filter threshold detect (0x0000 = lo, 0x0001 = hi)',
                'FIR1 threshold detect (0x0000 = lo, 0x0001 = hi)',
                'FIR2 threshold detect (0x0000 = lo, 0x0001 = hi)',
                'Feedback loop baseline compensation (x1000)',
                'Avg slope estimate over previous interval (x1000)',
                'Accumulated slope estimate over current interval (x1000)',
                'Slope update (x1000)']



#
# Read EDS Tool template file
# I assume that all defaults are filled in. What is not default is specified below
#
filename = 'eds_tool\\template.ini'
assert os.path.isfile(filename), "File  does not exist!"
print("Modify-ing ini file "+filename)
template_ini = open(filename, 'r')

filename = 'eds_tool\\525.ini'
new_ini = open(filename, 'w')






#
#
# HIER MOET JE WEZEN!!!
#
#
parameter_list = ['LpdA', 'LpdTau', 'Rc2', 'InputOffset', 'FIR1Threshold', 'FIR2Threshold',
                  'BcGainP', 'BcGainI', 'InputFile', 'TraceList', 'SizeLimit', 'EnergyResolution']

available_traces = [
    tracelist.index('Copy of input signal'),
    tracelist.index('HPD output (including datapath delay)'),
    tracelist.index('Bilateral filter output'),
    tracelist.index('FIR1 output (1st order)'),
    tracelist.index('FIR2 output (3rd order)'),
    tracelist.index('FIR1 threshold detect (0x0000 = lo, 0x0001 = hi)'),
    tracelist.index('FIR2 threshold detect (0x0000 = lo, 0x0001 = hi)'),
    tracelist.index('Slope update (x1000)')]

temp = str(available_traces[0])
n = 1
while n < len(available_traces):
    temp = temp + ',' + str(available_traces[n])
    n += 1

parameter_values = [-0.04, 150e-9, 5e-06, 0.9086, 0.0045, 0.003, 0, 0.02, '525.wav', temp, 1500000, 25e-6]
##parameter_values = [-0.3, 100e-9, 5e-06, 0.9002, 0, 0.02, 'lowcount.wav', temp, 1000000, 25e-6]




parameter_counter = 0
new_ini_data = []

for line in template_ini:
#    print line
    if parameter_counter < len(parameter_list):
        parameter = parameter_list[parameter_counter]

        # collect the traces from the ini file
        if line.find(parameter) > -1:
            if parameter == 'InputFile':
                print('Processing input file '+parameter_values[parameter_counter])

            new_ini.write(line[:-1] + str(parameter_values[parameter_counter]) + line[-1:])
            parameter_counter += 1
        else:
            new_ini.write(line)
    else:
        new_ini.write(line)

new_ini.close()


#
# Execute EDS Tool
#
os.chdir('eds_tool')
#os.remove('')
subprocess.call(['eds_toolV4.exe', '525.ini'])
os.chdir('..')



#
# Read wav file and convert # to Volt
#
filename_wavfile = 'eds_tool\\trace.wav'
(rate, data_value) = wavfile.read(filename_wavfile)
##data_value = (data_value.astype(np.float64)/2**15)  # ADC -32k.. +32k = -1V..1V
data_time = [x*0.025 for x in range(len(data_value[:,1]))]
print("Processing wav file "+filename_wavfile)


#
# prepare for plotting the result
#
trace_colour = ['k','r','g','b','c','m','y']

#
# locate the detected pulses in data_values, chop it up and average the chops
# this is based on the output of the FIR1 filter
#
selected_trace_number = [7,1]   # FIR1 output
pointer_available_traces = available_traces.index(selected_trace_number[0])
chopper = np.array(data_value[:,pointer_available_traces]).astype(np.bool)

# convolve FIR1 output with array of ones to suppress bounces in signal and to extend the width
n_taps = 601    # Okay for LPDA and Rc2
n_taps = 8001
chopper = sp.convolve(chopper, [1]*n_taps).astype(np.bool)
chopper = chopper[np.floor(n_taps/2):-np.floor(n_taps/2)]

# chop the HPD signal with the chopper signal
pointer_available_traces = available_traces.index(selected_trace_number[1])

# find the start of the chops in the chopped data
chops_pointer = np.concatenate(([0], chopper[1:]*-chopper[:-1]))
del chopper # clean up

# determine the number of chops, create a pointer and sort it
n_chops = sum(chops_pointer)
sum_chops = np.array([0]*n_taps)
chops_pointer = np.sort(chops_pointer * range(len(chops_pointer)))
chops_pointer = chops_pointer[-n_chops:]    # rhrow away the zeros in teh array


# First I check which step height occurs most often and determine the energy window from this.
step_height = []
nice_chop_counter = 0
temp = n_chops
while temp > 0:
    start = chops_pointer[len(chops_pointer)-temp]
    stop = start + n_taps

    left = np.average(data_value[start:start + np.floor(n_taps/4), pointer_available_traces])
    right = np.average(data_value[stop - np.floor(n_taps/4):stop, pointer_available_traces])

    if  (right-left) > -40000:      # I also consider negative slopes in my statistics
        step_height.append(right - left)

    temp -= 1

# Plot the histogram of step heights
fig = plt.figure("Histogram")
histogram_plot = fig.add_subplot(111)
histogram_plot.grid(True)
plt.xlabel('Step height [# of samples]')
plt.ylabel('Occurence')
plt.title('Distribution of step heights')
temp = get_current_fig_manager()
temp.window.SetPosition((0, 0))
temp.window.SetSize((300, 300))

(n, bins, temp) = histogram_plot.hist(step_height, bins=40)
del(temp)

# Determine the ranges of the energy window based on the bin with the largest # of occurences
pointer = 0
while n[pointer] < max(n):
    pointer += 1

energy_window = [bins[pointer], bins[pointer+1]]
print("The majority of steps has an energy range between "+str(bins[pointer])+" and "+str(bins[pointer+1]))


nice_chop_counter = 0
while n_chops > 0:
    start = chops_pointer[len(chops_pointer)-n_chops]
    stop = start + n_taps
    # skip resets
    left = np.average(data_value[start:start + np.floor(n_taps/4), pointer_available_traces])
    right = np.average(data_value[stop - np.floor(n_taps/4):stop, pointer_available_traces])

    # create an energy filter to suppress al steps that are to small/large
    timing_window = [-10,10]        # good for Rc2 and LPDA
##    timing_window = [-20,20]        # good for Offset
    d_threshold = 500

    if  (right-left) < energy_window[1] and (right-left) > energy_window[0]:
        # determine the 1st order derivative
        d_chopped_data = np.concatenate(([0], data_value[start+1:stop, pointer_available_traces] - data_value[start:stop-1, pointer_available_traces]))

        d_start = np.floor(n_taps/2) + 5 + timing_window[0] # this is where the slope starts
        d_stop = np.floor(n_taps/2) + timing_window[1]
        if max(d_chopped_data[d_start:d_stop]) > d_threshold:
            d_chopped_data_pointer = np.argmax(d_chopped_data[d_start:d_stop])-11
            sum_chops = sum_chops + data_value[start+d_chopped_data_pointer:stop+d_chopped_data_pointer, pointer_available_traces]
            nice_chop_counter +=1

    n_chops -= 1

del(d_chopped_data)

print("The total number of nice chops is "+str(nice_chop_counter))

# plot the chops
fig = plt.figure("Chops")
chops = fig.add_subplot(111)
chops.grid(True)
plt.xlabel('Time [us]')
plt.ylabel('I do not know')
plt.title('Average step')
temp = get_current_fig_manager()
temp.window.SetPosition((0, 300))
temp.window.SetSize((700, 750))


# Determine the slopes for LPDA, RC2 and Offset
# LPDA --> 200-500ns  after slope
# Rc2 --> 1-5us after slope
# Offset --> > 25us after slope
# system clock is 25ns
if nice_chop_counter > 0:   # if no chops are found, no slopes can be determined
    chops.plot(data_time[:n_taps], sum_chops, 'r.-', linewidth = 3, label = ("Average step"))

    # determine the middle of the slope
    d_sum_chops = np.concatenate(([0], sum_chops[1:]-sum_chops[:-1]))
    d_sum_chops_pointer = np.argmax(d_sum_chops)
    del d_sum_chops     # clean up

    #
    # plot slopes and selected region for LPDA
    #
    d_LPDA_range = np.array([200,500])/25  # time interval of interest [ns]

    # Determine the slope
    temp = sum_chops[d_sum_chops_pointer+d_LPDA_range[0]:d_sum_chops_pointer+d_LPDA_range[1]]
    d_LPDA = sp.stats.linregress(range(d_LPDA_range[1]-d_LPDA_range[0]),temp)
    print('d_LPDA = '+ str(d_LPDA[0]))

    # plot slopes and selected region for LPDA
    start = d_sum_chops_pointer + d_LPDA_range[0]
    stop = start - d_LPDA_range[0] + d_LPDA_range[1]
    d_range = np.array([-4,4]) * (d_LPDA_range[1]-d_LPDA_range[0])
    LPDA = np.array(range(d_range[0], d_range[1])) * d_LPDA[0]
    LPDA = LPDA + np.average(sum_chops[start:stop])
    chops.plot(data_time[(start+stop)/2 + d_range[0]:(start+stop)/2 + d_range[1]], LPDA, 'g', label = ("Slope LPDA setting"))
    chops.plot(data_time[start:stop], sum_chops[start:stop], 'g.-', linewidth = 5, label = ("LPDA slice"))

    del(LPDA)


    #
    # plot slopes and selected region for Rc2
    #
    d_Rc2_range = np.array([1000,5000])/25  # time interval of interest [ns]

    # Determine the slope
    temp = sum_chops[d_sum_chops_pointer+d_Rc2_range[0]:d_sum_chops_pointer+d_Rc2_range[1]]
    d_Rc2 = sp.stats.linregress(range(d_Rc2_range[1]-d_Rc2_range[0]),temp)
    print('d_Rc2 = '+ str(d_Rc2[0]))

    start = d_sum_chops_pointer + d_Rc2_range[0]
    stop = start - d_Rc2_range[0] + d_Rc2_range[1]
    d_range = np.array([-1,1]) * (d_Rc2_range[1]-d_Rc2_range[0])
    Rc2 = np.array(range(d_range[0], d_range[1])) * d_Rc2[0]
    Rc2 = Rc2 + np.average(sum_chops[start:stop])
    chops.plot(data_time[(start+stop)/2 + d_range[0]:(start+stop)/2 + d_range[1]], Rc2, 'k', label = ("Slope Rc2 setting"))
    chops.plot(data_time[start:stop], sum_chops[start:stop], 'k.-', linewidth = 5, label = ("Rc2 slice"))

    del(Rc2)



    #
    # plot slopes and selected region for InputOffset
    #
    d_IOffset_range = np.array([40000,80000])/25  # time interval of interest [ns]

    # Determine the slope
    temp = sum_chops[d_sum_chops_pointer+d_IOffset_range[0]:d_sum_chops_pointer+d_IOffset_range[1]]
    d_IOffset = sp.stats.linregress(range(d_IOffset_range[1]-d_IOffset_range[0]),temp)
    print('d_IOffset = '+ str(d_IOffset[0]))

    start = d_sum_chops_pointer + d_IOffset_range[0]
    stop = start - d_IOffset_range[0] + d_IOffset_range[1]
    d_range = np.array([-1,1]) * (d_IOffset_range[1]-d_IOffset_range[0])
    IOffset = np.array(range(d_range[0], d_range[1])) * d_IOffset[0]
    IOffset = IOffset + np.average(sum_chops[start:stop])
    chops.plot(data_time[(start+stop)/2 + d_range[0]:(start+stop)/2 + d_range[1]], IOffset, 'm', label = ("Slope InputOffset setting"))
    chops.plot(data_time[start:stop], sum_chops[start:stop], 'm.-', linewidth = 5, label = ("InputOffset slice"))

    # plot reference line with slope 0 for InputOffset, Rc2 and LPDA
    d_range = np.array([-2,1]) * (d_IOffset_range[1]-d_IOffset_range[0])
    ref_line = np.array(range(d_range[0], d_range[1])) * 0
    ref_line = ref_line + np.average(sum_chops[start:stop])
    chops.plot(data_time[(start+stop)/2 + d_range[0]:(start+stop)/2 + d_range[1]], ref_line, 'r--', label = ("Reference line"))

    del(IOffset)



    #
    # plot slopes and selected region for left plateau
    #
    d_left_plateau_range = np.array([-50000,-1000])/25  # time interval of interest [ns]

    # Determine the slope
    temp = sum_chops[d_sum_chops_pointer+d_left_plateau_range[0]:d_sum_chops_pointer+d_left_plateau_range[1]]
    d_left_plateau = sp.stats.linregress(range(d_left_plateau_range[1]-d_left_plateau_range[0]),temp)
    print('d_left_plateau = '+ str(d_left_plateau[0]))

    start = d_sum_chops_pointer + d_left_plateau_range[0]
    stop = start - d_left_plateau_range[0] + d_left_plateau_range[1]
    d_range = np.array([-1,1]) * (d_left_plateau_range[1]-d_left_plateau_range[0])
    left_plateau = np.array(range(d_range[0], d_range[1])) * d_left_plateau[0]
    left_plateau = left_plateau + np.average(sum_chops[start:stop])
    chops.plot(data_time[(start+stop)/2 + d_range[0]:(start+stop)/2 + d_range[1]], left_plateau, 'k', label = ("Slope left plateau"))
    chops.plot(data_time[start:stop], sum_chops[start:stop], 'k.-', linewidth = 5, label = ("Left plateau slice"))

    # plot reference line with slope 0 for left plateau
    d_range = np.array([-1,2]) * (d_left_plateau_range[1]-d_left_plateau_range[0])
    ref_line = np.array(range(d_range[0], d_range[1])) * 0
    ref_line = ref_line + np.average(sum_chops[start:stop])
    chops.plot(data_time[(start+stop)/2 + d_range[0]:(start+stop)/2 + d_range[1]], ref_line, 'r--', label = ("Reference line"))

    #
    # Legend
    #
    handles, labels = chops.get_legend_handles_labels()
    chops.legend(handles, labels)
    chops.legend(loc = 2)

del(left_plateau, sum_chops)






#
# This code plots sets of traces from the wave file
#
plot_range = np.array([0,160000])/25    # time in ns

fig1 = plt.figure()
HPD = fig1.add_subplot(111)
temp = get_current_fig_manager()
temp.window.SetPosition((700, 0))
temp.window.SetSize((700, 500))

HPD.set_xlabel('Time [us]')
HPD.set_ylabel('I dont know')
HPD.grid(True)

# plot HPD
selected_trace_number = [
    tracelist.index('Copy of input signal'),
    tracelist.index('HPD output (including datapath delay)'),
    tracelist.index('Bilateral filter output'),
    tracelist.index('FIR1 threshold detect (0x0000 = lo, 0x0001 = hi)')]

for x in range(len(selected_trace_number)):
    pointer_available_traces = available_traces.index(selected_trace_number[x])
    normalized_data = np.array(data_value[:,pointer_available_traces])
    normalized_data = normalized_data.astype(np.float64) / 2**15    #max(abs(normalized_data))
    HPD.plot(data_time, normalized_data, trace_colour[x], linewidth=1, label = (tracelist[selected_trace_number[x]]))

HPD.set_xlim(plot_range[0], plot_range[1])
HPD.set_ylim(-1, 1)

handles, labels = HPD.get_legend_handles_labels()
HPD.legend(handles, labels)

#
# Plot traces for FIR 1
#
fig2 = plt.figure("FIR1")
FIR1 = fig2.add_subplot(111)
temp = get_current_fig_manager()
temp.window.SetPosition((700, 500))
temp.window.SetSize((700, 500))

FIR1.set_xlabel('Time [us]')
FIR1.set_ylabel('Voltage [V]')
FIR1.grid(True)

selected_trace_number = [
    tracelist.index('HPD output (including datapath delay)')
    tracelist.index('Bilateral filter output'),
    tracelist.index('FIR1 output (1st order)'),
    tracelist.index('FIR1 threshold detect (0x0000 = lo, 0x0001 = hi)')]

FIR1_threshold = np.array(ones(len(data_time))) * parameter_values[parameter_list.index('FIR1Threshold')]
FIR1.plot(data_time, FIR1_threshold, 'r--', label = ("FIR1 threshold"))

for x in range(len(selected_trace_number)):
    pointer_available_traces = available_traces.index(selected_trace_number[x])
    normalized_data = np.array(data_value[:,pointer_available_traces])
    if x < len(selected_trace_number)-1:
        normalized_data = normalized_data.astype(np.float64) / 2**15    #max(abs(normalized_data))
    FIR1.plot(data_time, normalized_data, trace_colour[x], linewidth=1, label = (tracelist[selected_trace_number[x]]))

FIR1.set_xlim(plot_range[0], plot_range[1])
FIR1.set_ylim(-1, 1)

handles, labels = FIR1.get_legend_handles_labels()
FIR1.legend(handles, labels)


#
# Plot traces for FIR 2
#
fig3 = plt.figure("FIR2")
FIR2 = fig3.add_subplot(111)
temp = get_current_fig_manager()
temp.window.SetPosition((1000, 0))
temp.window.SetSize((700, 500))

FIR2.set_xlabel('Time [us]')
FIR2.set_ylabel('Voltage [V]')
FIR2.grid(True)
selected_trace_number = [
    tracelist.index('HPD output (including datapath delay)')
    tracelist.index('Bilateral filter output'),
    tracelist.index('FIR2 output (3rd order)'),
    tracelist.index('FIR2 threshold detect (0x0000 = lo, 0x0001 = hi)')]

FIR2_threshold = np.array(ones(len(data_time))) * parameter_values[parameter_list.index('FIR2Threshold')]
FIR2.plot(data_time, FIR2_threshold, 'r--', label = ("FIR2 threshold"))

for x in range(len(selected_trace_number)):
    pointer_available_traces = available_traces.index(selected_trace_number[x])
    normalized_data = np.array(data_value[:,pointer_available_traces])
    if x < len(selected_trace_number)-1:
        normalized_data = normalized_data.astype(np.float64) / 2**15    #max(abs(normalized_data))
    FIR2.plot(data_time, normalized_data, trace_colour[x], linewidth=1, label = (tracelist[selected_trace_number[x]]))

FIR2.set_xlim(plot_range[0], plot_range[1])
FIR2.set_ylim(-1, 1)

handles, labels = FIR2.get_legend_handles_labels()
FIR2.legend(handles, labels)



plt.show()
