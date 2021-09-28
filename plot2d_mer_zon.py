#!/usr/bin/env python

import glob
import numpy as n
import matplotlib.pyplot as plt
import h5py
import cfi_dac as cfi

name="all"
#name="summer_mer_zon_50km_900s"

#name="winter_mer_zon_100km_900s"

lag_name="Meridional lag (km)"
lag_name2="Zonal lag (km)"
fl=glob.glob("mpi/%s/*.h5"%(name))
fl.sort()
print(fl)
h=h5py.File(fl[0],"r")
s_x=n.copy(h["s_x"][()])
s_y=n.copy(h["s_y"][()])
dtau=n.copy(h["dtau"][()])
ds_x=n.copy(h["ds_x"][()])

acfs=n.zeros([6,len(s_x),len(s_y)])
ws=n.zeros([6,len(s_x),len(s_y)])
h.close()
n_avg=0.0
for f in fl:
    print(f)
    h=h5py.File(f,"r")
    if "acf" in h.keys():
        if n.sum(n.isnan(h["acf"][()])) == 0:
            # print(h["err"][()])
            w=1.0#/h["err"][()]
            xi=h["xi"][()]
            yi=h["yi"][()]
            acfs[:,xi,yi]+=w*h["acf"][()]
            ws[:,xi,yi]+=w
    h.close()
acfs=acfs/ws
print(s_x)
dx=s_x[1]-s_x[0]
xx=n.concatenate((s_x-dx/2,[n.max(s_x)+dx/2.0]))
yy=n.concatenate((s_y-dx/2,[n.max(s_y)+dx/2.0]))
print(xx)


plt.figure(figsize=(8,8))
plt.subplot(321)

#acfs[0,4,0]=0.0

plt.pcolormesh(xx,yy,acfs[0,:,:],vmin=-0.2e3,vmax=1.5e3)
plt.title(cfi.cf_names[0])
plt.xlabel(lag_name2)
plt.ylabel(lag_name)
plt.colorbar()
plt.subplot(322)

plt.pcolormesh(xx,yy,acfs[1,:,:],vmin=-0.2e3,vmax=1.5e3)
plt.title(cfi.cf_names[1])
plt.xlabel(lag_name2)
plt.ylabel(lag_name)
plt.colorbar()
plt.subplot(323)

plt.pcolormesh(xx,yy,acfs[2,:,:],vmin=-100,vmax=100)
plt.title(cfi.cf_names[2])
plt.xlabel(lag_name2)
plt.ylabel(lag_name)
plt.colorbar()

plt.subplot(324)
plt.pcolormesh(xx,yy,acfs[3,:,:],vmin=-200,vmax=200)
plt.title(cfi.cf_names[3])
plt.xlabel(lag_name2)
plt.ylabel(lag_name)
plt.colorbar()

plt.subplot(325)
plt.pcolormesh(xx,yy,acfs[2,:,:],vmin=-100,vmax=100)
plt.title(cfi.cf_names[4])
plt.xlabel(lag_name2)
plt.ylabel(lag_name)
plt.colorbar()

plt.subplot(326)
plt.pcolormesh(xx,yy,acfs[2,:,:],vmin=-100,vmax=100)
plt.title(cfi.cf_names[5])
plt.xlabel(lag_name2)
plt.ylabel(lag_name)
plt.colorbar()
plt.tight_layout()
plt.show()
