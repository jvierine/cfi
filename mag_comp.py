#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as n

import supermag_read as smr
import mmaria_read as mmr
import cfi_config as c
import scipy.interpolate as sio

def remove_baseline(zon):
    """
    Remove all low frequency (more than 6 hour period) waves
    """
    bidx=n.where(n.isnan(zon))[0]
    for i in bidx:
        for j in range(24):
            if not n.isnan(zon[i+j]):
                zon[i]=zon[i+j]
                break
    bix=n.where(n.isnan(zon))[0]
    print(bix)
    Z=n.fft.fft(zon)
    f=n.fft.fftfreq(len(Z),d=1800.0)
    low_idx=n.where(n.abs(f) < 1.0/(6.0*3600.0))
    Z[low_idx]=0.0
    zonb=n.real(n.fft.ifft(Z))
    return(zonb)

def lpf_mag(mag,N=60):
    return(n.convolve(n.abs(mag)**2.0,n.repeat(1.0/N,N),mode="same"))

d=smr.supermag_data(fname="tro_2018_2019.csv",plot=False)
b=d.get_bounds()
print(b)
# get all data from "TRO"
m=d.get_meas(b[0],b[0]+365*24*3600,station="TRO")


# directory with all mmaria network data
md=mmr.mmaria_data(c.data_directory)
# what is the data bounds (first and last time stamp)
print(md.get_bounds())

# read all meteor radar data between these two timestamps
d=md.read_data(b[0],b[0]+365*24*3600,read_all_detections=False)

zon=d["v"][0,:,:]
mer=d["v"][1,:,:]

#zonb=remove_baseline(zon[:,20])
zonb=remove_baseline(zon[:,25])
merb=remove_baseline(mer[:,25])
wken=zonb**2.0+merb**2.0

plt.plot(d["times"],zonb)
print(d["times"][1]-d["times"][0])
plt.show()

mag=n.copy(n.sqrt(m["neu"][:,0]**2.0+m["neu"][:,1]**2.0))
mag=lpf_mag(mag)


zf=sio.interp1d(d["times"],wken,bounds_error=False)

H,x,y=n.histogram2d(n.log10(zf(m["times"])),n.log10(mag),bins=100,range=[[-2,5],[-2,7]])
plt.pcolormesh(x,y,n.transpose(n.log10(H+1.0)))
plt.xlabel("Mesospheric wind kinetic energy ($m^2/s^2$)")
plt.ylabel("Magnetometer proxy for Joule heating ($(nT)^2$)")
plt.show()

#plt.pcolormesh(d["times"],d["rgs"]/1e3,n.transpose(d["v"][0,:,:]),vmin=-100,vmax=100)
#plt.show()
