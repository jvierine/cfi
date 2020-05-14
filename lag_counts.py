#!/usr/bin/env python
#
# Estimate the number of horizontal and temporal lags
# 
import h5py
import numpy as n
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpi4py import MPI

import cfi_config as c
import mmaria_read as mr
import geoid_const as gc
import stuffr

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

latdeg2km=gc.latdeg2km
londeg2km=gc.londeg2km

#works for mmaria_data
md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()

t0=stuffr.date2unix(2019,6,1,0,0,0)
t1=stuffr.date2unix(2019,6,15,0,0,0)

# about two weeks of data in June 2019
d=md.read_data(t0,t1)

dcos_thresh=0.8
lats=d["lats"]
lons=d["lons"]
heights=d["heights"]            #hights in meters
heights=heights/1000
dcoss=d["dcos"]
t=d["t"]

dcos2=n.sqrt(dcoss[:,0]**2.0+dcoss[:,1]**2.0)
ok_idx=n.where(dcos2 < dcos_thresh)[0]

lats=lats[ok_idx]
lons=lons[ok_idx]
heights=heights[ok_idx]
t=t[ok_idx]

print("total number of measurements %d"%(len(t)))

# define a histogram
n_m=len(lats)
rgs=n.arange(80,105)
dh=1.0
h_bins=100

# base 10 log horizontal distance
mer_dist_bins=n.linspace(-2,3.2,num=h_bins+1)
zon_dist_bins=n.linspace(-2,3.2,num=h_bins+1)
hor_dist_bins=n.linspace(-2,3.2,num=h_bins+1)
ver_dist_bins=n.linspace(-50,50,num=h_bins+1)

# base 10 log temporal distance
tdist_bins=n.linspace(-2,6.2,num=h_bins+1)

# zonal, meridional, horizontal, and vertical distances
HZ=n.zeros([h_bins,h_bins])
HM=n.zeros([h_bins,h_bins])
HH=n.zeros([h_bins,h_bins])
HV=n.zeros([h_bins,h_bins])

# parallel processing.
# processor number=comm.rank
# number if processes=comm.size
for ri in range(comm.rank,len(rgs),comm.size):
    HZ[:,:]=0.0
    HM[:,:]=0.0
    HH[:,:]=0.0
    HV[:,:]=0.0
    rg=rgs[ri]

    zon_dists=[]
    mer_dists=[]
    hor_dists=[]    
    ver_dists=[]
    t_dists=[]

    max_zon_dist=0
    max_mer_dist=0
    max_hor_dist=0
    
    idx0=n.where( (heights >= rg) & (heights < (rg+dh)) )[0]
    print("number of measurements between %1.2f and %1.2f km = %d"%(rg,rg+dh,len(idx0)))

    lats0 = n.copy(lats[idx0])
    lons0 = n.copy(lons[idx0])
    t0s   = n.copy(t[idx0])
    h0s   = n.copy(heights[idx0])
    
    for mi in range(len(lats0)):
        if mi % 100 == 0:
            print("rank %d %d/%d"%(comm.rank,mi,len(idx0)))
        lat0=lats0[mi]
        lon0=lons0[mi]
        t0=t0s[mi]
        h0=h0s[mi]        


        # horizontal distances
        hor_dists=n.sqrt((n.abs(lats0[mi:len(lats0)]-lat0)*latdeg2km)**2.0 + (n.abs(lons0[mi:len(lats0)]-lon0)*londeg2km)**2.0 )
        mer_dists=n.abs(lats0[mi:len(lats0)]-lat0)*latdeg2km
        zon_dists=n.abs(lons0[mi:len(lats0)]-lon0)*londeg2km
        ver_dists=h0s[mi:len(h0s)]-h0

        mmd=n.max(mer_dists)
        if mmd > max_mer_dist:
            max_mer_dist=mmd
        mzd=n.max(zon_dists)
        if mzd > max_zon_dist:
            max_zon_dist=mzd
        mhd=n.max(hor_dists)
        if mhd > max_mer_dist:
            max_hor_dist=mhd
            
        temporal_dists = n.abs(t[mi:len(lats0)]-t0)
        HZ0,zxe,zye=n.histogram2d(n.log10(zon_dists+0.01),n.log10(temporal_dists+0.01),bins=[hor_dist_bins,tdist_bins])
        HM0,mxe,mye=n.histogram2d(n.log10(mer_dists+0.01),n.log10(temporal_dists+0.01),bins=[hor_dist_bins,tdist_bins])
        HH0,hxe,hye=n.histogram2d(n.log10(hor_dists+0.01),n.log10(temporal_dists+0.01),bins=[hor_dist_bins,tdist_bins])
        HZ+=HZ0
        HM+=HM0
        HH+=HH0

    plt.pcolormesh(mxe,mye,n.transpose(n.log10(HH+1)))
    cb=plt.colorbar()
    cb.set_label("$\log_{10}$ (counts)")
    plt.xlabel("Horizontal distance $\log_{10}$ (km)")
    plt.ylabel("Temporal distance $\log_{10}$ (s)")
    plt.title("Distribution of horizontal and temporal lags, h=%1.0f km\nMax distance=%1.0f km"%(rg,max_hor_dist))
    plt.tight_layout()
    plt.savefig("./figs/fig_hor_%03d.png"%(rg))
    plt.close()
    plt.clf()

    plt.pcolormesh(zxe,zye,n.transpose(n.log10(HZ+1)))
    cb=plt.colorbar()
    cb.set_label("$\log_{10}$ (counts)")
    plt.xlabel("Zonal distance $\log_{10}$ (km)")
    plt.ylabel("Temporal distance $\log_{10}$ (s)")
    plt.title("Distribution of zonal and temporal lags, h=%1.0f km\nMax distance=%1.0f km"%(rg,max_zon_dist))
    plt.tight_layout()
    plt.savefig("./figs/fig_zon_%03d.png"%(rg))
    plt.close()
    plt.clf()

    plt.pcolormesh(mxe,mye,n.transpose(n.log10(HM+1)))
    cb=plt.colorbar()
    cb.set_label("$\log_{10}$ (counts)")
    plt.xlabel("Meridional distance $\log_{10}$ (km)")
    plt.ylabel("Temporal distance $\log_{10}$ (s)")
    plt.title("Distribution of meridional and temporal lags, h=%1.0f km\nMax distance=%1.0f km"%(rg,max_mer_dist))
    plt.tight_layout()
    plt.savefig("./figs/fig_mer_%03d.png"%(rg))
    plt.close()
    plt.clf()

    
        
        
    ho=h5py.File("hist_%1.2f.h5"%(rg),"w")
    ho["time_diff_log10_s"]=tdist_bins
    ho["horizontal_diff_log10_km"]=hor_dist_bins
    ho["H_hor"]=HH
    ho["H_zon"]=HZ
    ho["H_mer"]=HM
    ho["tbins"]=tdist_bins
    ho["hx"]=hxe
    ho["hy"]=hye    
    ho["zx"]=zxe
    ho["zy"]=zye    
    ho["mx"]=mxe
    ho["my"]=mye    
    ho.close()

