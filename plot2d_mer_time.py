#!/usr/bin/env python

import glob
import numpy as n
import matplotlib.pyplot as plt
import h5py
name="summer_2d_mer_time"
fl=glob.glob("mpi/%s/*.h5"%(name))
fl.sort()
h=h5py.File(fl[0],"r")
acfs=n.copy(h["acfs"].value)
errs=n.copy(h["errs"].value)
tau=n.copy(h["tau"].value)
s_y=n.copy(h["s_y"].value)
acfs[:,:,:]=0.0
errs[:,:,:]=0.0
ws=n.copy(acfs)
h.close()
n_avg=0.0
for f in fl:
    print(f)
    h=h5py.File(f,"r")
    if "acfs" in h.keys():
        w=1.0/h["errs"].value
        acfs+=w*h["acfs"].value
        ws+=w
        n_avg+=1.0
    h.close()
acfs=acfs/ws
plt.pcolormesh(tau,s_y,acfs[0,:,:],vmin=0,vmax=1000)
plt.title("$C_{uu}$")
plt.xlabel("Time (s)")
plt.ylabel("Meridional distance (km)")
plt.colorbar()
plt.show()

plt.pcolormesh(tau,s_y,acfs[1,:,:],vmin=0,vmax=1000)
plt.title("$C_{vv}$")
plt.xlabel("Time (s)")
plt.ylabel("Meridional distance (km)")

plt.colorbar()
plt.show()

plt.pcolormesh(tau,s_y,acfs[2,:,:],vmin=-200,vmax=200)
plt.title("$C_{uv}$")
plt.xlabel("Time (s)")
plt.ylabel("Meridional distance (km)")

plt.colorbar()
plt.show()
