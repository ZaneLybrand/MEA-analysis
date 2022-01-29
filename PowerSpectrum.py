'''
Max peak for FFT                                                                            Lybrand Lab 2021
The function of this code is to find the average power spectral density (PSD) across multiple wells of an MEA recording.

Below you can input the file name of the data in .csv format located in the parent python file. When executed, this code
will return:
1. A six panel figure that shows the average of 2 groups of raw signal, the average FFT of the defined wells, and
spectrogram of the raw signals.
2. A second figure displays the individual FFT of all activated wells.
3. A .csv will be created in the parent python folder that reports all peaks in the individual FFT above determined
threshold (adjustable below).
4. A list of maximum peak FFT values from activated wells will be displayed in the command prompt. This should
correspond to the peak values returned in the .csv.

For questions and comments, contact Zane Lybrand (zlybrand@twu.edu)
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal
plt.style.use('default')

data = pd.read_csv("250kpa_5000hz_compare.csv")  #Change to file for analysis
t = data["time"]
A1 = data["A1"]
A2 = data["A2"]
A3 = data["A3"]
A4 = data["A4"]
A5 = data["A5"]
A6 = data["A6"]
B1 = data["B1"]
B2 = data["B2"]
B3 = data["B3"]
B4 = data["B4"]
B5 = data["B5"]
B6 = data["B6"]
C1 = data["C1"]
C2 = data["C2"]
C3 = data["C3"]
C4 = data["C4"]
C5 = data["C5"]
C6 = data["C6"]
D1 = data["D1"]
D2 = data["D2"]
D3 = data["D3"]
D4 = data["D4"]
D5 = data["D5"]
D6 = data["D6"]

### Power spectral density (PSD)
dt = t[1]-t[0]
N= A1.shape[0]
T = N * dt

### Calculate Power spectral density (PSD) for all channels

xfA1 = np.fft.fft(A1 - A1.mean())                   # Compute FFT of A1
A1_s = 2 * dt ** 2 / T * (xfA1 * np.conj(xfA1))     # Compute Spectrum
A1_s = A1_s[:int(len(A1) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQA1 = 1 / dt / 2                                   # Nyquist frequency
faxis = np.arange(0, fNQA1, df)                      # Frequency axis

xfA2 = np.fft.fft(A2 - A2.mean())                   # Compute FFT of A2
A2_s = 2 * dt ** 2 / T * (xfA2 * np.conj(xfA2))     # Compute Spectrum
A2_s = A2_s[:int(len(A2) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQA2 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQA2, df)                     # Frequency axis

xfA3 = np.fft.fft(A3 - A3.mean())                   # Compute FFT of A3
A3_s = 2 * dt ** 2 / T * (xfA3 * np.conj(xfA3))     # Compute Spectrum
A3_s = A3_s[:int(len(A3) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQA3 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQA3, df)                     # Frequency axis

xfA4 = np.fft.fft(A4 - A4.mean())                   # Compute FFT of A4
A4_s = 2 * dt ** 2 / T * (xfA4 * np.conj(xfA4))     # Compute Spectrum
A4_s = A4_s[:int(len(A4) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQA4 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQA4, df)                     # Frequency axis

xfA5 = np.fft.fft(A5 - A5.mean())                   # Compute FFT of A5
A5_s = 2 * dt ** 2 / T * (xfA5 * np.conj(xfA5))     # Compute Spectrum
A5_s = A5_s[:int(len(A5) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQA5 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQA5, df)                     # Frequency axis

xfA6 = np.fft.fft(A6 - A6.mean())                   # Compute FFT of A6
A6_s = 2 * dt ** 2 / T * (xfA6 * np.conj(xfA6))     # Compute Spectrum
A6_s = A6_s[:int(len(A6) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQA6 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQA6, df)                     # Frequency axis

xfB1 = np.fft.fft(B1 - B1.mean())                   # Compute FFT of B1
B1_s = 2 * dt ** 2 / T * (xfB1 * np.conj(xfB1))     # Compute Spectrum
B1_s = B1_s[:int(len(B1) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQB1 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQB1, df)                     # Frequency axis

xfB2 = np.fft.fft(B2 - B2.mean())                   # Compute FFT of B2
B2_s = 2 * dt ** 2 / T * (xfB2 * np.conj(xfB2))     # Compute Spectrum
B2_s = B2_s[:int(len(B2) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQB2 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQB2, df)                     # Frequency axis

xfB3 = np.fft.fft(B3 - B3.mean())                   # Compute FFT of B3
B3_s = 2 * dt ** 2 / T * (xfB3 * np.conj(xfB3))     # Compute Spectrum
B3_s = B3_s[:int(len(B3) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQB3 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQB3, df)                     # Frequency axis

xfB4 = np.fft.fft(B4 - B4.mean())                   # Compute FFT of B4
B4_s = 2 * dt ** 2 / T * (xfB4 * np.conj(xfB4))     # Compute Spectrum
B4_s = B4_s[:int(len(B4) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQB4 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQB4, df)                     # Frequency axis

xfB5 = np.fft.fft(B5 - B5.mean())                   # Compute FFT of B5
B5_s = 2 * dt ** 2 / T * (xfB5 * np.conj(xfB5))     # Compute Spectrum
B5_s = B5_s[:int(len(B5) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQB5 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQB5, df)                     # Frequency axis

xfB6 = np.fft.fft(B6 - B6.mean())                   # Compute FFT of B6
B6_s = 2 * dt ** 2 / T * (xfB6 * np.conj(xfB6))     # Compute Spectrum
B6_s = B6_s[:int(len(B6) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQB6 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQB6, df)

xfC1 = np.fft.fft(C1 - C1.mean())                   # Compute FFT of C1
C1_s = 2 * dt ** 2 / T * (xfC1 * np.conj(xfC1))     # Compute Spectrum
C1_s = C1_s[:int(len(C1) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQC1 = 1 / dt / 2                                   # Nyquist frequency
faxis = np.arange(0, fNQC1, df)                      # Frequency axis

xfC2 = np.fft.fft(C2 - C2.mean())                   # Compute FFT of C2
C2_s = 2 * dt ** 2 / T * (xfC2 * np.conj(xfC2))     # Compute Spectrum
C2_s = C2_s[:int(len(C2) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQC2 = 1 / dt / 2                                   # Nyquist frequency
faxis = np.arange(0, fNQC2, df)                      # Frequency axis

xfC3 = np.fft.fft(C3 - C3.mean())                   # Compute FFT of C3
C3_s = 2 * dt ** 2 / T * (xfC3 * np.conj(xfC3))     # Compute Spectrum
C3_s = C3_s[:int(len(C3) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQC3 = 1 / dt / 2                                   # Nyquist frequency
faxis = np.arange(0, fNQC3, df)                      # Frequency axis

xfC4 = np.fft.fft(C4 - C4.mean())                   # Compute FFT of C4
C4_s = 2 * dt ** 2 / T * (xfC4 * np.conj(xfC4))     # Compute Spectrum of
C4_s = C4_s[:int(len(C4) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQC4 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQC4, df)                     # Frequency axis

xfC5 = np.fft.fft(C5 - C5.mean())                   # Compute FFT of C5
C5_s = 2 * dt ** 2 / T * (xfC5 * np.conj(xfC5))     # Compute Spectrum of
C5_s = C5_s[:int(len(C5) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQC5 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQC5, df)                     # Frequency axis

xfC6 = np.fft.fft(C6 - C6.mean())                   # Compute FFT of C6
C6_s = 2 * dt ** 2 / T * (xfC6 * np.conj(xfC6))     # Compute Spectrum
C6_s = C6_s[:int(len(C6) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQC6 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQC6, df)

xfD1 = np.fft.fft(D1 - D1.mean())                   # Compute FFT of D1
D1_s = 2 * dt ** 2 / T * (xfD1 * np.conj(xfD1))     # Compute Spectrum
D1_s = D1_s[:int(len(D1) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQD1 = 1 / dt / 2                                   # Nyquist frequency
faxis = np.arange(0, fNQD1, df)                      # Frequency axis

xfD2 = np.fft.fft(D2 - D2.mean())                   # Compute FFT of D2
D2_s = 2 * dt ** 2 / T * (xfD2 * np.conj(xfD2))     # Compute Spectrum
D2_s = D2_s[:int(len(D2) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQD2 = 1 / dt / 2                                   # Nyquist frequency
faxis = np.arange(0, fNQD2, df)                      # Frequency axis

xfD3 = np.fft.fft(D3 - D3.mean())                   # Compute FFT of D3
D3_s = 2 * dt ** 2 / T * (xfD3 * np.conj(xfD3))     # Compute Spectrum
D3_s = D3_s[:int(len(D3) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQD3 = 1 / dt / 2                                   # Nyquist frequency
faxis = np.arange(0, fNQD3, df)                      # Frequency axis

xfD4 = np.fft.fft(D4 - D4.mean())                   # Compute FFT of D4
D4_s = 2 * dt ** 2 / T * (xfD4 * np.conj(xfD4))     # Compute Spectrum
D4_s = D4_s[:int(len(D4) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQD4 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQD4, df)                     # Frequency axis

xfD5 = np.fft.fft(D5 - D5.mean())                   # Compute FFT of D5
D5_s = 2 * dt ** 2 / T * (xfD5 * np.conj(xfD5))     # Compute Spectrum
D5_s = D5_s[:int(len(D5) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQD5 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQD5, df)                     # Frequency axis

xfD6 = np.fft.fft(D6 - D6.mean())                   # Compute FFT of D6
D6_s = 2 * dt ** 2 / T * (xfD6 * np.conj(xfD6))     # Compute Spectrum
D6_s = D6_s[:int(len(D6) / 2)]                      # Remove negative frequencies
df = 1 / T.max()                                    # Determine frequency resolution
fNQD6 = 1 / dt / 2                                  # Nyquist frequency
faxis = np.arange(0, fNQD6, df)


''' In this section you can define your groups for FFT by changing which WellID_s is in each group. 
Currently it can calculate 2 different groups (Healthy vs Epileptic) with the output as an average FFT across the group.
You will be provided the peak values and the maximal peak for each group'''
## Compile average FFT across groups
g1FFT = np.array([A1_s, A2_s, A3_s, A4_s, A5_s, A6_s, B1_s, B2_s, B3_s, B4_s, B5_s, B6_s])      #change this line for group 1 (include _s)

g2FFT = np.array([C1_s, C2_s, C3_s, C4_s, C5_s, C6_s, D1_s, D2_s, D3_s, D4_s, D5_s, D6_s])       #change this line for group 2

g1_avgFFT = np.average(g1FFT,0)
g2_avgFFT = np.average(g2FFT,0)

## Define average array for raw data
group1 = np.array([A1,A2,A3,A4,A5,A6,B1,B2,B3,B4,B5,B6])                                   #change this line for group 1
group2 = np.array([C1,C2,C3,C4,C5,C6,D1,D2,D3,D4,D5,D6])                                    #change this line for group 2

g1_avg = np.average(group1,0)
g2_avg = np.average(group2,0)

### Find Peaks of FFT
th = .01
s1peak = signal.find_peaks(np.real(g1_avgFFT), prominence=th)
s2peak = signal.find_peaks(np.real(g2_avgFFT), prominence=th)

A1peak = signal.find_peaks(np.real(A1_s), prominence=th)
A2peak = signal.find_peaks(np.real(A2_s), prominence=th)
A3peak = signal.find_peaks(np.real(A3_s), prominence=th)
A4peak = signal.find_peaks(np.real(A4_s), prominence=th)
A5peak = signal.find_peaks(np.real(A5_s), prominence=th)
A6peak = signal.find_peaks(np.real(A6_s), prominence=th)
B1peak = signal.find_peaks(np.real(B1_s), prominence=th)
B2peak = signal.find_peaks(np.real(B2_s), prominence=th)
B3peak = signal.find_peaks(np.real(B3_s), prominence=th)
B4peak = signal.find_peaks(np.real(B4_s), prominence=th)
B5peak = signal.find_peaks(np.real(B5_s), prominence=th)
B6peak = signal.find_peaks(np.real(B6_s), prominence=th)
C1peak = signal.find_peaks(np.real(C1_s), prominence=th)
C2peak = signal.find_peaks(np.real(C2_s), prominence=th)
C3peak = signal.find_peaks(np.real(C3_s), prominence=th)
C4peak = signal.find_peaks(np.real(C4_s), prominence=th)
C5peak = signal.find_peaks(np.real(C5_s), prominence=th)
C6peak = signal.find_peaks(np.real(C6_s), prominence=th)
D1peak = signal.find_peaks(np.real(D1_s), prominence=th)
D2peak = signal.find_peaks(np.real(D2_s), prominence=th)
D3peak = signal.find_peaks(np.real(D3_s), prominence=th)
D4peak = signal.find_peaks(np.real(D4_s), prominence=th)
D5peak = signal.find_peaks(np.real(D5_s), prominence=th)
D6peak = signal.find_peaks(np.real(D6_s), prominence=th)

# g1peak = signal.find_peaks(np.real(g1FFT), prominence=.05)
# g2peak = signal.find_peaks(np.real(g2FFT), prominence=.05)

s1heighty = s1peak[1]['prominences']
# s1heightx = s1peak[0]

s2heighty = s2peak[1]['prominences']
# s2heightx = s2peak[0]

fftmax1 = np.amax(s1heighty)
fftmax2 = np.amax(s2heighty)

#find Max Peak of FFT for each channel
A1max = np.amax(A1peak[1]['prominences'], initial=0)
A2max = np.amax(A2peak[1]['prominences'], initial=0)
A3max = np.amax(A3peak[1]['prominences'], initial=0)
A4max = np.amax(A4peak[1]['prominences'], initial=0)
A5max = np.amax(A5peak[1]['prominences'], initial=0)
A6max = np.amax(A6peak[1]['prominences'], initial=0)
B1max = np.amax(B1peak[1]['prominences'], initial=0)
B2max = np.amax(B2peak[1]['prominences'], initial=0)
B3max = np.amax(B3peak[1]['prominences'], initial=0)
B4max = np.amax(B4peak[1]['prominences'], initial=0)
B5max = np.amax(B5peak[1]['prominences'], initial=0)
B6max = np.amax(B6peak[1]['prominences'], initial=0)
C1max = np.amax(C1peak[1]['prominences'], initial=0)
C2max = np.amax(C2peak[1]['prominences'], initial=0)
C3max = np.amax(C3peak[1]['prominences'], initial=0)
C4max = np.amax(C4peak[1]['prominences'], initial=0)
C5max = np.amax(C5peak[1]['prominences'], initial=0)
C6max = np.amax(C6peak[1]['prominences'], initial=0)
D1max = np.amax(D1peak[1]['prominences'], initial=0)
D2max = np.amax(D2peak[1]['prominences'], initial=0)
D3max = np.amax(D3peak[1]['prominences'], initial=0)
D4max = np.amax(D4peak[1]['prominences'], initial=0)
D5max = np.amax(D5peak[1]['prominences'], initial=0)
D6max = np.amax(D6peak[1]['prominences'], initial=0)


#Convert Peak to array
A1array = np.array(A1peak[1]['prominences'])
A2array = np.array(A2peak[1]['prominences'])
A3array = np.array(A3peak[1]['prominences'])
A4array = np.array(A4peak[1]['prominences'])
A5array = np.array(A5peak[1]['prominences'])
A6array = np.array(A6peak[1]['prominences'])
B1array = np.array(B1peak[1]['prominences'])
B2array = np.array(B2peak[1]['prominences'])
B3array = np.array(B3peak[1]['prominences'])
B4array = np.array(B4peak[1]['prominences'])
B5array = np.array(B5peak[1]['prominences'])
B6array = np.array(B6peak[1]['prominences'])
C1array = np.array(C1peak[1]['prominences'])
C2array = np.array(C2peak[1]['prominences'])
C3array = np.array(C3peak[1]['prominences'])
C4array = np.array(C4peak[1]['prominences'])
C5array = np.array(C5peak[1]['prominences'])
C6array = np.array(C6peak[1]['prominences'])
D1array = np.array(D1peak[1]['prominences'])
D2array = np.array(D2peak[1]['prominences'])
D3array = np.array(D3peak[1]['prominences'])
D4array = np.array(D4peak[1]['prominences'])
D5array = np.array(D5peak[1]['prominences'])
D6array = np.array(D6peak[1]['prominences'])


""" This Portion creates an array of the max peak from individual FFT and saves as .csv"""
MaxArray = np.array([A1max, A2max, A3max, A4max,A5max,A6max, B1max, B2max, B3max,B4max,B5max,B6max, C1max, C2max, C3max, C4max, C5max, C6max, D1max, D2max, D3max, D4max, D5max, D6max ])

dfA1 = pd.DataFrame({'A1':A1array})
dfA2 = pd.DataFrame({'A2':A2array})
dfA3 = pd.DataFrame({'A3':A3array})
dfA4 = pd.DataFrame({'A4':A4array})
dfA5 = pd.DataFrame({'A5':A5array})
dfA6 = pd.DataFrame({'A6':A6array})

dfB1 = pd.DataFrame({'B1':B1array})
dfB2 = pd.DataFrame({'B2':B2array})
dfB3 = pd.DataFrame({'B3':B3array})
dfB4 = pd.DataFrame({'B4':B4array})
dfB5 = pd.DataFrame({'B5':B5array})
dfB6 = pd.DataFrame({'B6':B6array})

dfC1 = pd.DataFrame({'C1':C1array})
dfC2 = pd.DataFrame({'C2':C2array})
dfC3 = pd.DataFrame({'C3':C3array})
dfC4 = pd.DataFrame({'C4':C4array})
dfC5 = pd.DataFrame({'C5':C5array})
dfC6 = pd.DataFrame({'C6':C6array})

dfD1 = pd.DataFrame({'D1':D1array})
dfD2 = pd.DataFrame({'D2':D2array})
dfD3 = pd.DataFrame({'D3':D3array})
dfD4 = pd.DataFrame({'D4':D4array})
dfD5 = pd.DataFrame({'D5':D5array})
dfD6 = pd.DataFrame({'D6':D6array})
df = pd.concat([dfA1,dfA2, dfA3,dfA4,dfA5,dfA6, dfB1,dfB2, dfB3,dfB4,dfB5,dfB6, dfC1,dfC2, dfC3, dfC4, dfC5, dfC6, dfD1,dfD2, dfD3, dfD4, dfD5, dfD6], ignore_index=True, axis=1)    #Adds NaN to empty objects to make all columns equal length

df.to_csv("Compiled Array.csv", index=True, header=True)


print(" A1 max peak = ", A1max)
print(" A2 max peak = ", A2max)
print(" A3 max peak = ", A3max)
print(" A4 max peak = ", A4max)
print(" A5 max peak = ", A5max)
print(" A6 max peak = ", A6max)
print(" B1 max peak = ", B1max)
print(" B2 max peak = ", B2max)
print(" B3 max peak = ", B3max)
print(" B4 max peak = ", B4max)
print(" B5 max peak = ", B5max)
print(" B6 max peak = ", B6max)
print(" C1 max peak = ", C1max)
print(" C2 max peak = ", C2max)
print(" C3 max peak = ", C3max)
print(" C4 max peak = ", C4max)
print(" C5 max peak = ", C5max)
print(" C6 max peak = ", C6max)
print(" D1 max peak = ", D1max)
print(" D2 max peak = ", D2max)
print(" D3 max peak = ", D3max)
print(" D4 max peak = ", D4max)
print(" D5 max peak = ", D5max)
print(" D6 max peak = ", D6max)


""" This section sets up the plot for each in a 2X3 configuration"""

### Building Spectrogram
Fs = 1 / dt                                      # Sampling frequency
interval = int(Fs)                               # interval size
overlap = int(Fs * 0.95)                         # overlap interval

### Build plots
plt.subplot(231)
plt.plot(t, g1_avg)
plt.title('1hr LFP')
plt.ylim([-10,10])                              # change this for y-axis limit on Group 1 avg plot
plt.xlabel('Time (s)')
plt.ylabel(r'Voltage($\mu V$)')

plt.subplot(232)
plt.plot(faxis, np.real(g1_avgFFT))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.ylim([0,5])
plt.title('Group 1 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(233)
f, t1, Sxx1 = signal.spectrogram(
    g1_avg,
    fs=Fs,
    nperseg=interval,
    noverlap=overlap)
plt.pcolormesh(t1, f, 10 * np.log10(Sxx1),
           cmap='jet', shading='auto')
plt.colorbar()
plt.ylim([0,100])
plt.xlabel('Time[s]')
plt.ylabel('Frequency [Hz]')

plt.subplot(234)
plt.plot(t, g2_avg, color='#ff7f0e')
plt.ylim([-10,10])
plt.title('24hr LFP')
plt.xlabel('Time (s)')
plt.ylabel(r'Voltage($\mu V$)')

plt.subplot(235)
plt.plot(faxis, np.real(g2_avgFFT))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.ylim([0,5])
plt.title('Group 2 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(236)
f, t2, Sxx2 = signal.spectrogram(
    g2_avg,
    fs=Fs,
    nperseg=interval,
    noverlap=overlap)
plt.pcolormesh(t2, f, 10 * np.log10(Sxx2),
           cmap='jet', shading='auto')
plt.colorbar()
plt.ylim([0,100])
plt.xlabel('Time[s]')
plt.ylabel('Frequency [Hz]')


plt.tight_layout()
plt.show()
""" This portion will plot a separate figure to overlay FFT"""
plt.figure()
plt.plot(faxis, np.real(g1_avgFFT), linewidth=4)                     #Plot spectrum vs freq
plt.plot(faxis, np.real(g2_avgFFT), linewidth=4)
plt.xlim([0,30])
plt.ylim([0,30])
plt.title('250kPa 5000Hz PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')


plt.tight_layout()
plt.show()
"""This section will plot all individual FFT signals"""
plt.figure()
plt.subplot(231)
plt.plot(faxis, np.real(A1_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('A1 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(232)
plt.plot(faxis, np.real(A2_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('A2 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(233)
plt.plot(faxis, np.real(A3_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('A3 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(234)
plt.plot(faxis, np.real(A4_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('A4 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(235)
plt.plot(faxis, np.real(A5_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('A5 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(236)
plt.plot(faxis, np.real(A6_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('A6 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.figure()
plt.subplot(231)
plt.plot(faxis, np.real(B1_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('B1 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(232)
plt.plot(faxis, np.real(B2_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('B2 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(233)
plt.plot(faxis, np.real(B3_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('B3 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(234)
plt.plot(faxis, np.real(B4_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('B4 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(235)
plt.plot(faxis, np.real(B5_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('B5 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(236)
plt.plot(faxis, np.real(B6_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('B6 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.figure()
plt.subplot(231)
plt.plot(faxis, np.real(C1_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('C1 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(232)
plt.plot(faxis, np.real(C2_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('C2 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(233)
plt.plot(faxis, np.real(C3_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('C3 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(234)
plt.plot(faxis, np.real(C4_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('C4 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(235)
plt.plot(faxis, np.real(C5_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('C5 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(236)
plt.plot(faxis, np.real(C6_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('C6 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.figure()
plt.subplot(231)
plt.plot(faxis, np.real(D1_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('D1 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(232)
plt.plot(faxis, np.real(D2_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('D2 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(233)
plt.plot(faxis, np.real(D3_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('D3 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(234)
plt.plot(faxis, np.real(D4_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('D4 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(235)
plt.plot(faxis, np.real(D5_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('D5 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.subplot(236)
plt.plot(faxis, np.real(D6_s))                     #Plot spectrum vs freq
plt.xlim([0,100])
plt.title('D6 PSD')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Power [$\mu V^2$/Hz]')

plt.tight_layout()
plt.show()
