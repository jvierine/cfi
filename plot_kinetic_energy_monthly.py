#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import h5py
import glob

alt=90

run_name="small_50km_ke"

fl=glob.glob("mpi/00.00_%s/*.h5"%(run_name))
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


dl=glob.glob("mpi/??.??_%s"%(run_name))
dl.sort()
print(dl)
n_days=len(dl)

uu=n.zeros([acfs.shape[1],n_days])
vv=n.zeros([acfs.shape[1],n_days])


tods=n.concatenate((tods,[24.0]))
for mi in range(n_days):
    dirn=dl[mi]
    fl=glob.glob("%s/*.h5"%(dirn))
    fl.sort()
    acfs[:,:,:]=0.0
    n_avg=0.0
    ws=n.zeros(errs.shape)
    
    for f in fl:
        print(f)
        h=h5py.File(f,"r")
        gidx=n.where(n.isnan(h["acfs"].value)!=True)
        w=n.zeros(errs.shape)    
        #w[gidx]=1.0/h["errs"].value[gidx]
        #ws[gidx]+=w[gidx]
        #acfs[gidx]+=w[gidx]*h["acfs"].value[gidx]

        w=1.0#n.zeros(errs.shape)    
        ws+=w
        acfs+=w*h["acfs"].value
        h.close()
        n_avg+=1.0
    acfs=acfs/ws
    uu[:,mi]=acfs[0,:,hi]
    vv[:,mi]=acfs[1,:,hi]
#    plt.pcolormesh(tods,heights,n.log10(n.transpose(acfs[1,:,:]+acfs[0,:,:])),vmin=0,vmax=4)
 #   plt.colorbar()
  #  plt.show()


months=n.linspace(0,12,num=n_days+1)

plt.figure(figsize=(12,6))
plt.subplot(131)
plt.title("Meridional")
plt.pcolormesh(tods,months,n.transpose(uu),vmin=0,vmax=6000)
plt.colorbar()
plt.xlabel("Hour of day (UTC)")
plt.ylabel("Month of year")

plt.subplot(132)
plt.title("Zonal")
plt.pcolormesh(tods,months,n.transpose(vv),vmin=0,vmax=6000)
plt.colorbar()
plt.xlabel("Hour of day (UTC)")
plt.ylabel("Month of year")

plt.subplot(133)
plt.title("Total")
plt.pcolormesh(tods,months,n.log10(n.transpose(vv+uu)),vmin=2,vmax=4)
plt.colorbar()
plt.xlabel("Hour of day (UTC)")
plt.ylabel("Month of year")

plt.show()



plt.pcolormesh(n.log10( vv+uu),vmin=2,vmax=4)
plt.show()
