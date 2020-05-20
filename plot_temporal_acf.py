#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import glob
import scipy.signal as ss
import h5py
import cfi_dac as cfi

name="summer_tacf_koki"
#name="winter_tacf_koki"
fl=glob.glob("mpi/%s/*.h5"%(name))
fl.sort()
h=h5py.File(fl[0],"r")
tau=n.copy(h["tau"].value)
acfs=n.copy(h["acfs"].value)
acfs[:,:]=0.0
ws=n.copy(h["acfs"].value)
ws[:,:]=0.0
print(ws)
print(acfs)

h.close()
n_avg=1.0
for f in fl:
    h=h5py.File(f,"r")
    w=1/(h["errs"].value)

    wacf=w*n.copy(h["acfs"].value)
    # print(n.sum(wacf))
    # print(n.sum(n.isnan(wacf)))
    
  #  print(n.sum(w))
   # print(n.sum(n.isnan(w))    )
    if (n.sum(n.isnan(wacf))) == 0:
        acfs+=wacf
        ws+=w


    n_avg+=1.0
    if False:
        plt.plot(tau,w[:,0])
        plt.show()
        
        plt.plot(tau,h["acfs"].value[:,0])
        plt.plot(tau,h["acfs"].value[:,1])
        plt.show()

    h.close()
print(ws)
acfs=acfs/ws
wf=ss.hann(len(tau)*2)[len(tau):(2*len(tau))]
plt.subplot(121)
plt.plot(tau/(24*3600),acfs[:,0])
plt.plot(tau/(24*3600),acfs[:,1])
plt.xlabel("Time lag (days)")
plt.ylabel("Correlation function (m$^2$/s$^2$)")

plt.subplot(122)
wf=1.0

Suu=n.fft.hfft(wf*acfs[:,0])
Svv=n.fft.hfft(wf*acfs[:,1])

dt=n.diff(tau)[0]
f=n.fft.fftfreq(len(Suu),d=dt)
#print(Suu)
Suu[0]=0
fi=n.argmax(Suu)
print(fi)
print(Suu[fi])
f0=n.abs(f[fi])
print(f0)
#b=Suu[fi]/f0**(a)

plt.loglog(n.real(n.fft.fftshift(f)),n.fft.fftshift(Suu))
plt.loglog(n.real(n.fft.fftshift(f)),n.fft.fftshift(Svv))

b=Suu[fi]/f0**(-5/3.0)
print(b)
plt.loglog(n.abs(n.fft.fftshift(f)),0.1*b*n.abs(n.fft.fftshift(f))**(-5/3.0),label="-5/3",alpha=0.3)
b=Suu[fi]/f0**(-2)
plt.loglog(n.abs(n.fft.fftshift(f)),0.1*b*n.abs(n.fft.fftshift(f))**(-2),label="-2",alpha=0.3)
plt.legend()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power spectral density (m$^2$/s$^2$/Hz)")
plt.axvline(1.0/(24*3600.0),color="gray",linestyle="--")
plt.axvline(1.0/(12*3600.0),color="gray",linestyle="--")
plt.axvline(1.0/(8*3600.0),color="gray",linestyle="--")




plt.show()
    
    
