#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import h5py
import glob

import msis_pglow as mp

#name="large_ke"
#name="small_50km_ke"
name="test_ke"
# winter months
winter=True
if winter:
    m0=12
    title="Winter"
    fl0=glob.glob("mpi/10.??_%s/*.h5"%(name))
    fl1=glob.glob("mpi/11.??_%s/*.h5"%(name))
    fl2=glob.glob("mpi/00.??_%s/*.h5"%(name))
else:
    # summer months
    title="Summer"
    m0=6
    fl0=glob.glob("mpi/05.??_%s/*.h5"%(name))
    fl1=glob.glob("mpi/06.??_%s/*.h5"%(name))
    fl2=glob.glob("mpi/07.??_%s/*.h5"%(name))

fl=fl0+fl1+fl2
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
mass=n.zeros([len(tods),len(heights)])
for ti,t in enumerate(tods):
    for hi,h in enumerate(heights):
        mass[ti,hi]=mp.get_atmospheric_mass(2010,m0,hour=int(n.floor(t)),alt=h)

tods=n.concatenate((tods,[23.99]))
print(tods)

plt.figure(figsize=(8,8))
plt.subplot(221)
plt.pcolormesh(tods,heights,n.transpose(acfs[0,:,:]),vmin=0,vmax=8000)
plt.colorbar()
plt.xlabel("Time of day (UTC hour)")
plt.ylabel("Height (km)")
plt.title("%s $G_{uu}$"%(title))
plt.subplot(222)

plt.pcolormesh(tods,heights,n.transpose(acfs[1,:,:]),vmin=0,vmax=8000)
plt.xlabel("Time of day (UTC hour)")
plt.ylabel("Height (km)")
plt.title("$G_{vv}$")
plt.colorbar()
plt.subplot(223)

plt.pcolormesh(tods,heights,n.transpose(acfs[2,:,:]),vmin=0,vmax=1000)
plt.xlabel("Time of day (UTC hour)")
plt.ylabel("Height (km)")
plt.title("$G_{ww}$")
plt.colorbar()
plt.subplot(224)

plt.pcolormesh(tods,heights,n.log10(0.5*n.transpose(mass)*n.transpose(acfs[1,:,:]+acfs[0,:,:]+acfs[2,:,:])),vmin=-5,vmax=-1)
plt.xlabel("Time of day (UTC hour)")
plt.ylabel("Height (km)")
plt.title("$\\frac{1}{2}\\rho(G_{uu}+G_{vv}+G_{ww})$ (J/m$^2$)" )
plt.colorbar()
plt.tight_layout()
plt.show()

    
