''' Spike visualization                                      Lybrand Lab 2021
The function of this code is to plot the average spikes and overlay standard deviation
across 2 different MEA channels
Below you can input the file name of the data in .csv format located in the parent python file. When executed this code
will:
1. Compute an average trace for all input events.
2. Plot a two panel figure comparing shape and amplitude of events with standard deviation of all input events

For input format (.csv) see companion file SpikeExample1 and SpikeExample2

For questions and comments, contact Zane Lybrand (zlybrand@twu.edu)'''
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal
plt.style.use('default')

data = pd.read_csv("SpikeExample1.csv")      #Access File 1
data = np.array(data)
t = data[:,0]                                       #define time axis
n = len(data[0])                                    #determine number of spikes
spikes = data[:,1:(n-1)]                            #define spikes from data
avg_spike = np.average(spikes,axis=1)               #compute average spikes
std_spike = np.std(spikes, axis=1)                  #compute stdev of spikes
upstd = avg_spike + std_spike
dwnstd = avg_spike - std_spike


data2 = pd.read_csv("SpikeExample1.csv")  #Access File 2
data2 = np.array(data2)
t2 = data2[:,0]
n2 = len(data2[0])
spikes2 = data2[:,1:(n2-1)]
avg_spike2 = np.average(spikes2,axis=1)
std_spike2 = np.std(spikes2, axis=1)
upstd2 = avg_spike2 + std_spike2
dwnstd2= avg_spike2 - std_spike2

plt.subplot(211)
plt.plot(t, avg_spike, label="1hr")
plt.fill_between(t, dwnstd, upstd, facecolor="lightsteelblue")
plt.title('250kPa at 5000Hz Population Spike')
# plt.ylim([-10,10])                       #change for specifying y-axis limit
plt.xlabel('Time (ms)')
plt.ylabel(r'Voltage($\mu V$)')
# plt.legend(frameon=False)                #change to add legend

plt.subplot(212)
plt.plot(t2, avg_spike2, color='darkorange', label="24hr")
plt.fill_between(t2, dwnstd2, upstd2, facecolor="moccasin")
# plt.title("250kPa at 3000Hz Population Spike")
# plt.ylim([-10,10])                        #change for y-axis
# plt.xlim([-1,10])                         #change for x-axis
plt.xlabel('Time (ms)')
plt.ylabel(r'Voltage($\mu V$)')
# plt.legend(frameon=False)

plt.tight_layout()
plt.show()


