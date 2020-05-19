#!/usr/bin/env python

import numpy as n
import glob
import matplotlib.pyplot as plt
import h5py
import scipy.signal as ss

h0=95
fl=glob.glob("/mnt/data/juha/mmaria_norway/mean_wind/*.h5")
fl.sort()
h=h5py.File(fl[0],"r")
v=n.copy(h["v"].value)
heights=n.copy(h["rgs"].value)
times=n.copy(h["times"].value)
n_t=len(times)
time_h=24.0*n.arange(n_t)/n_t

hi=n.argmin(n.abs(heights-h0))

n_days=len(fl)

U=n.zeros([n_t,n_days])
V=n.zeros([n_t,n_days])

h.close()
u=n.array([])
v=n.array([])
dudx=n.array([])
dudy=n.array([])
dvdx=n.array([])
dvdy=n.array([])
vort=n.array([])
hdiv=n.array([])
q=[]
for di in range(n_days):
    h=h5py.File(fl[di],"r")
    U[:,di]=h["v"].value[0,:,hi]
    u=n.concatenate((u,h["v"].value[0,:,hi]))
    v=n.concatenate((v,h["v"].value[1,:,hi]))
    dudy=n.concatenate((dudy,h["v"].value[2,:,hi]))
    dvdy=n.concatenate((dvdy,h["v"].value[3,:,hi]))
    dudx=n.concatenate((dudx,h["v"].value[4,:,hi]))
    dvdx=n.concatenate((dvdx,h["v"].value[5,:,hi]))

    dudyp=h["v"].value[2,:,hi]
    dvdxp=h["v"].value[5,:,hi]
    vort=n.concatenate((vort, -dudyp+dvdxp))
    
    dudxp=h["v"].value[4,:,hi]
    dvdyp=h["v"].value[3,:,hi]
    hdiv=n.concatenate((hdiv,dudxp+dvdyp))
    

    V[:,di]=h["v"].value[1,:,hi]

u[n.isnan(u)]=0.0
v[n.isnan(v)]=0.0
dudx[n.isnan(dudx)]=0.0
dudy[n.isnan(dudy)]=0.0
dvdx[n.isnan(dvdx)]=0.0
dvdy[n.isnan(dvdy)]=0.0
vort[n.isnan(vort)]=0.0
hdiv[n.isnan(hdiv)]=0.0

#vort=n.real(n.fft.ifft(n.fft.fft(n.repeat(1.0/4.0,4.0),len(vort))*n.fft.fft(vort)))
#hdiv=n.real(n.fft.ifft(n.fft.fft(n.repeat(1.0/4.0,4.0),len(vort))*n.fft.fft(hdiv)))

f=n.fft.fftshift(n.fft.fftfreq(len(u),d=900.0))

def filterlp(x,L):
    return(n.real(n.fft.ifft(n.fft.fft(n.repeat(1.0/float(L),L),len(x))*n.fft.fft(x))))
    
plt.subplot(321)

plt.axvline(1.0/(24*3600),color="gray")
plt.axvline(1.0/(12*3600),color="gray")
plt.axvline(1.0/(8*3600),color="gray")
plt.axvline(1.0/(27.32*24*3600),color="gray")
plt.axvline(1.0/(365*24*3600),color="gray")
plt.loglog(n.abs(f),filterlp(n.abs(n.fft.fftshift(n.fft.fft(u)))**2.0,10),color="red")
#plt.loglog(n.abs(f),n.abs(f)**(-3))
#plt.loglog(n.abs(f),n.abs(f)**(-5/3.0))
plt.subplot(322)
plt.axvline(1.0/(24*3600),color="gray")
plt.axvline(1.0/(12*3600),color="gray")
plt.axvline(1.0/(8*3600),color="gray")
plt.axvline(1.0/(27.32*24*3600),color="gray")
plt.axvline(1.0/(365*24*3600),color="gray")
plt.loglog(n.abs(f),filterlp(n.abs(n.fft.fftshift(n.fft.fft(v)))**2.0,10),color="red")
plt.subplot(323)
plt.axvline(1.0/(24*3600),color="gray")
plt.axvline(1.0/(12*3600),color="gray")
plt.axvline(1.0/(8*3600),color="gray")
plt.axvline(1.0/(6*3600),color="gray")
plt.axvline(1.0/(27.32*24*3600),color="gray")
plt.axvline(1.0/(365*24*3600),color="gray")
plt.loglog(n.abs(f),filterlp(n.abs(n.fft.fftshift(n.fft.fft(vort)))**2.0,10),color="red")
plt.subplot(324)
plt.axvline(1.0/(24*3600),color="gray")
plt.axvline(1.0/(12*3600),color="gray")
plt.axvline(1.0/(8*3600),color="gray")
plt.axvline(1.0/(6*3600),color="gray")
plt.axvline(1.0/(27.32*24*3600),color="gray")
plt.axvline(1.0/(365*24*3600),color="gray")
plt.loglog(n.abs(f),filterlp(n.abs(n.fft.fftshift(n.fft.fft(hdiv)))**2.0,10),color="red")
plt.subplot(325)
plt.axvline(1.0/(24*3600),color="gray")
plt.axvline(1.0/(12*3600),color="gray")
plt.axvline(1.0/(8*3600),color="gray")
plt.axvline(1.0/(27.32*24*3600),color="gray")
plt.axvline(1.0/(365*24*3600),color="gray")
plt.loglog(n.abs(f),n.abs(n.fft.fftshift(n.fft.fft(dvdx)))**2.0,color="red")
plt.subplot(326)
plt.axvline(1.0/(24*3600),color="gray")
plt.axvline(1.0/(12*3600),color="gray")
plt.axvline(1.0/(8*3600),color="gray")
plt.axvline(1.0/(27.32*24*3600),color="gray")
plt.axvline(1.0/(365*24*3600),color="gray")
plt.loglog(n.abs(f),n.abs(n.fft.fftshift(n.fft.fft(dvdy)))**2.0,color="red")
plt.show()
    
U[n.isnan(U)]=0.0
V[n.isnan(V)]=0.0

plt.subplot(311)
plt.title("Zonal h=%1.0f km"%(h0))
plt.pcolormesh(U,vmin=-150,vmax=150)
plt.ylabel("Time of day (h)")
plt.xlabel("Day")

plt.subplot(312)
plt.title("Meridional")
plt.pcolormesh(V,vmin=-150,vmax=150)
plt.ylabel("Time of day (h)")
plt.xlabel("Day")

plt.subplot(313)
KE=V**2.0+U**2.0
for i in range(KE.shape[0]):
    KE[i,:]=n.roll(n.real(n.fft.ifft(n.fft.fft(n.repeat(1.0/31.0,31),n_days)*n.fft.fft(KE[i,:]))),-15)
plt.title("Monthly mean kinetic energy")
plt.pcolormesh(n.log10(KE),vmin=0,vmax=4)
plt.ylabel("Time of day (h)")
plt.xlabel("Day")
plt.show()


