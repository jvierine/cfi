#!/usr/bin/env python

import glob
import numpy as n
import matplotlib.pyplot as plt
import h5py
import cfi_dac as cfi

#name="summer_2d_zon_time_3"
#lag_name="Zonal distance"

name="summer_2d_mer_time_12"
lag_name="Meridional distance"

#name="winter_2d_mer_time_12_900s_100km"
#lag_name="Meridional distance"

#name="winter_2d_mer_time_12_900s_50km"
#lag_name="Meridional distance"

#name="winter_2d_zon_time_12_900s_50km"
#lag_name="Meridional distance"

fl=glob.glob("mpi/%s/*.h5"%(name))
fl.sort()

h=h5py.File(fl[0],"r")
y_lags=n.copy(h["y_lags"].value)
t_lags=n.copy(h["t_lags"].value)
n_y=len(y_lags)
n_t=len(t_lags)
acfs=n.zeros([6,n_y,n_t])
ws=n.zeros([6,n_y,n_t])
h.close()

for f in fl:
    print(f)
    h=h5py.File(f,"r")
    if "acfs" in h.keys():
        print(h["acfs"].value)
        if n.sum(n.isnan(h["acfs"].value)) == 0:
            yi=h["yi"].value
            ti=h["ti"].value
            w=1.0/h["errs"].value
            acfs[:,yi,ti]+=w*h["acfs"].value
            ws[:,yi,ti]+=w
    h.close()
acfs=acfs/ws

xtau=n.linspace(n.min(t_lags),n.max(t_lags),num=len(t_lags)+1)
yy=n.linspace(n.min(y_lags),n.max(y_lags),num=len(y_lags)+1)

plt.figure(figsize=(8,8))
plt.subplot(321)


plt.pcolormesh(xtau/3600.0,yy,acfs[0,:,:],vmin=-0.2e3,vmax=1.5e3)
plt.title(cfi.cf_names[0])
plt.xlabel("Time (h)")
plt.ylabel("%s (km)"%(lag_name))
plt.colorbar()
plt.subplot(322)

plt.pcolormesh(xtau/3600.0,yy,acfs[1,:,:],vmin=-0.5e3,vmax=1.5e3)
plt.title(cfi.cf_names[1])
plt.xlabel("Time (h)")
plt.ylabel("%s (km)"%(lag_name))
plt.colorbar()
plt.subplot(323)

plt.pcolormesh(xtau/3600.0,yy,acfs[2,:,:],vmin=-100,vmax=100)
plt.title(cfi.cf_names[2])
plt.xlabel("Time (h)")
plt.ylabel("%s (km)"%(lag_name))
plt.colorbar()

plt.subplot(324)
plt.pcolormesh(xtau/3600.0,yy,acfs[3,:,:],vmin=-200,vmax=200)
plt.title(cfi.cf_names[3])
plt.xlabel("Time (h)")
plt.ylabel("%s (km)"%(lag_name))
plt.colorbar()

plt.subplot(325)
plt.pcolormesh(xtau/3600.0,yy,acfs[2,:,:],vmin=-100,vmax=100)
plt.title(cfi.cf_names[4])
plt.xlabel("Time (h)")
plt.ylabel("%s (km)"%(lag_name))
plt.colorbar()

plt.subplot(326)
plt.pcolormesh(xtau/3600.0,yy,acfs[2,:,:],vmin=-100,vmax=100)
plt.title(cfi.cf_names[5])
plt.xlabel("Time (h)")
plt.ylabel("%s (km)"%(lag_name))
plt.colorbar()
plt.tight_layout()
plt.show()
