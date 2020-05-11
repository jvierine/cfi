#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import h5py
import glob


fl=glob.glob("mpi/winter_small_ke/*.h5")
fl.sort()

h=h5py.File(fl[0],"r")
acfs=n.copy(h["acfs"])
errs=n.copy(h["errs"].value)
tods=n.copy(h["tods"].value)
heights=n.copy(h["heights"].value)
ds_h=h["ds_h"].value
acfs[:,:,:]=0.0
h.close()
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
tods=n.concatenate((tods,[24.0]))
plt.pcolormesh(tods,heights,n.transpose(acfs[0,:,:]),vmin=0,vmax=8000)
plt.colorbar()
plt.show()
plt.pcolormesh(tods,heights,n.transpose(acfs[1,:,:]),vmin=0,vmax=8000)
plt.colorbar()
plt.show()

plt.pcolormesh(tods,heights,n.log10(n.transpose(acfs[1,:,:]+acfs[0,:,:])),vmin=0,vmax=4)
plt.colorbar()
plt.show()

    
