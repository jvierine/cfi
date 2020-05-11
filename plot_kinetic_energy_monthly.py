#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import h5py
import glob

alt=90.0

run_name="large_ke"

fl=glob.glob("mpi/00_%s/*.h5"%(run_name))
fl.sort()

h=h5py.File(fl[0],"r")
acfs=n.copy(h["acfs"])
errs=n.copy(h["errs"].value)
tods=n.copy(h["tods"].value)
heights=n.copy(h["heights"].value)
ds_h=h["ds_h"].value
acfs[:,:,:]=0.0
h.close()
#print(acfs.shape)

hi=n.argmin(n.abs(alt-heights))

uu=n.zeros([acfs.shape[1],12])
vv=n.zeros([acfs.shape[1],12])
tods=n.concatenate((tods,[24.0]))
for mi in range(12):
    fl=glob.glob("mpi/%02d_%s/*.h5"%(mi,run_name))
    fl.sort()
    acfs[:,:,:]=0.0
    n_avg=0.0
    ws=n.zeros(errs.shape)
    for f in fl:
        print(f)
        h=h5py.File(f,"r")
        w=1.0/h["errs"].value
        if n.sum(n.isnan(h["acfs"].value)) > 0:
            ws+=w
            acfs+=w*h["acfs"].value
        h.close()
        n_avg+=1.0
    acfs=acfs/ws
    uu[:,mi]=acfs[0,:,hi]
    vv[:,mi]=acfs[1,:,hi]
    plt.pcolormesh(tods,heights,n.log10(n.transpose(acfs[1,:,:]+acfs[0,:,:])),vmin=0,vmax=4)
    plt.colorbar()
    plt.show()

    

plt.pcolormesh(tods,n.arange(13),n.transpose(uu),vmin=0,vmax=8000)
plt.colorbar()
plt.show()
plt.pcolormesh(tods,n.arange(13),n.transpose(vv),vmin=0,vmax=8000)
plt.colorbar()
plt.show()
plt.pcolormesh(tods,n.arange(13),n.log10(n.transpose(vv+uu)),vmin=0,vmax=4)
plt.colorbar()
plt.show()
