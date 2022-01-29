''' WELCOME USER                                                             Lybrand Lab 2021
This code is used to calculate the phase angle difference between two chosen electrodes extrapolated from
multielectrode arrays. The data are sampled from .csv files found in the Python Analysis folder.

For each file to analyze, (1) change the file name and (2) select the channels of interest.

When executed this code will return 2 figures and a value for the phase angle synchronization.
 
First figure will display:
1) A plot with the raw signals from the two selected channels
2) A plot overlaying both channels with the Morlet Wavelet transformation
3) A plot overlaying the phase angle time series for both channels

The second figure will display each individual phase angles (grey) on a radial plot. Overlayed will be the average
phase angle (red).

For example of formatting data file (.csv) please see sample data file. 

For any questions of comments please contact Zane Lybrand (zlybrand@twu.edu) 
'''
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use('default')

data = pd.read_csv("2020118_High_5kHz_001_D6_combo.csv")  #Change file name for analysis
t = data["time"]
ch21 = data["ch21"]
ch31 = data["ch31"]
ch12 = data["ch12"]
ch22 = data["ch22"]
ch32 = data["ch32"]
ch42 = data["ch42"]
ch13 = data["ch13"]
ch23 = data["ch23"]
ch33 = data["ch33"]
ch43 = data["ch43"]
ch24 = data["ch24"]
ch34 = data["ch34"]


### Select Channels to compare                          #Change to channels of interest (e.g. ch12, ch21, etc.)
s1 = ch22
s2 = ch21

## Set up FFT and convolution

# Create complex Morlet Wavelet
srate = len(t)
center_freq = 5                                         # filter frequency
time = np.linspace(-1, 1, srate)

# Create complex sine wave
sine_wave = np.exp((1j*2) * np.pi * center_freq * time)

# Create Gaussian window
s = 7 / (2 * np.pi * center_freq)
gaus_win = np.exp((-np.power(time,2) / (2 * np.power(s,2))))

# Create Morlet Wavelet
cmw = np.multiply(sine_wave, gaus_win)
half_wavN = (len(time)-1)/2

#FFT of wavelet
n_wavelet = len(cmw)
n_data    = len(t)
n_conv    = n_wavelet + (n_data-1)

#FFT of wavelet
waveletX = np.fft.fft(cmw)
max_waveletX = np.amax(waveletX)
waveletXf = np.divide(waveletX, max_waveletX)

# compute Hz for plotting
hz_min = 0.0
hz_max = np.floor((n_conv/2)+1)
hz_half = np.multiply(srate, 0.5)
hz = np.linspace(hz_min, hz_max, len(t))

# plot wavelet
# plt.figure(1)
# plt.subplot(211)
# plt.plot(time, cmw.real)
# plt.title('Morlet Wavelet')
# plt.xlim(-1,1)
#
# plt.subplot(212)
# plt.plot(hz, waveletX)
# plt.title('Power Spectrum of Morlet Wavelet')
# plt.xlim(0, 50)

# plt.tight_layout()
# plt.show()

## Calculate phase angle for each channel
# analytic signal of channel 1
fft_data = np.fft.fft(s1)                                   #fft of channel 1
ch1_conv = np.multiply(waveletX, fft_data)                  #convolution of channel 1
an_s1 = np.fft.ifft(ch1_conv)                               #inverse fft to return analog signal

# collect real and phase data
phase_data1 = np.angle(an_s1)
real_data1 = np.real(an_s1)

# analytic signal of channel 2
fft_data = np.fft.fft(s2)                                   #fft of channel 2
ch2_conv = np.multiply(waveletX, fft_data)                  #convolution of channel 2
an_s2 = np.fft.ifft(ch2_conv)                               #inverse fft to return analog signal

# # collect real and phase data
phase_data2 = np.angle(an_s2)
real_data2 = np.real(an_s2)

## SET UP FIGURE AND PLOT
plt.figure(2)
plt.subplot(223)
plt.plot(t, real_data1)
plt.plot(t, real_data2)
# plt.xlim(0, 2)
plt.title('Filtered signal')
plt.xlabel('Time(s)')
plt.ylabel(r'Voltage($\mu V$)')

plt.subplot(224)
plt.plot(t, phase_data1)
plt.plot(t, phase_data2)
# plt.xlim(0, 2)
plt.title('Phase angle time series')
plt.xlabel('Time(s)')
plt.ylabel(r'Phase angle (radian)')

plt.subplot(221)
plt.plot(t, s1)
plt.title('Channel 1' )
plt.xlabel('Time(s)')
plt.ylabel(r'Voltage($\mu V$)')
plt.ylim(-20, 20)

plt.subplot(222)
plt.plot(t, s2, '#F97306')
plt.title('Channel 2')
plt.xlabel('Time(s)')
plt.ylabel(r'Voltage($\mu V$)')
plt.ylim(-20, 20)

### CALCULATE PHASE ANGLE DIFFERENCE
phase_angle_differences = phase_data2-phase_data1               #phase angle difference
euler_phase_differences = np.exp(1j*phase_angle_differences)    #euler representation of angles
mean_complex_vector = np.mean(euler_phase_differences)          #mean vector (in complex space)
phase_synchronization = np.absolute(mean_complex_vector)        #length of mean vector (M from Me^ik)
print('Synchronization between ', phase_synchronization)

#Plot phase angles in polar space
plt.figure(3)

theta = phase_angle_differences
r = np.tile([0, 1], 25000)

ax = plt.subplot(111, projection='polar')
ax.plot(theta, r, linewidth=0.5, color='grey')
ax.plot([0, np.angle(mean_complex_vector)], [0, phase_synchronization], linewidth=3, color='firebrick')

ax.set_rmax(1)
ax.set_rticks([0.25, 0.5, 0.75, 1.0])  # Less radial ticks
ax.set_rlabel_position(-22.5)          # Move radial labels away from plotted line
ax.grid(True)
ax.set_title("Phase angle difference", va='bottom')

plt.tight_layout()
plt.show()
