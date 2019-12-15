#!/usr/bin/env python
import h5py
import numpy as n
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

pairs=[]
latdeg2km=111.321
londeg2km=65.122785

h=h5py.File("res/simone_nov2018_multilink_juha_30min_1000m.h5","r")
lats=h["lats"].value
lons=h["lons"].value
heights=h["heights"].value
t=h["t"].value
print("total number of measurements %d"%(len(t)))

# define a histogram
n_m=len(lats)
rgs=n.arange(80,100)
dh=1.0

h_bins=100
# base 10 log horizontal distance
hdist_bins=n.linspace(-1,3.2,num=h_bins+1)
# base 10 log temporal distance
tdist_bins=n.linspace(-1,6.2,num=h_bins+1)
H=n.zeros([h_bins,h_bins])

#for ri,rg in enumerate(rgs):
# parallel processing.
# processor number=comm.rank
# number if processes=comm.size
for ri in range(comm.rank,len(rgs),comm.size):
    H[:,:]=0.0
    rg=rgs[ri]

    hor_dists=[]
    ver_dists=[]
    t_dists=[]

    idx0=n.where( (heights >= rg) & (heights < (rg+dh)) )[0]
    print("number of measurements between %1.2f and %1.2f km = %d"%(rg,rg+dh,len(idx0)))

    lats0 = n.copy(lats[idx0])
    lons0 = n.copy(lons[idx0])
    t0s   = n.copy(t[idx0])
    
    for mi in range(len(lats0)):
        if mi % 100 == 0:
            print("rank %d %d/%d"%(comm.rank,mi,len(idx0)))
        lat0=lats0[mi]
        lon0=lons0[mi]
        t0=t0s[mi]            

        # horizontal distances
        hor_dists=n.sqrt((n.abs(lats0[mi:len(lats0)]-lat0)*latdeg2km)**2.0 + (n.abs(lons0[mi:len(lats0)]-lon0)*londeg2km)**2.0 )
#        print(n.max(hor_dists))
 #       print(n.log10(n.max(hor_dists)))
        # temporal distances
#        plt.plot(hor_dists,".")
 #       plt.show()
        temporal_dists = n.abs(t[mi:len(lats0)]-t0)
        H0,xe,ye=n.histogram2d(n.log10(hor_dists+0.1),n.log10(temporal_dists+0.1),bins=[hdist_bins,tdist_bins])
        H+=H0
        #print(n.max(xe))
        #print(n.max(ye))        
    plt.pcolormesh(xe,ye,n.transpose(n.log10(H+1)))
    cb=plt.colorbar()
    cb.set_label("$\log_{10}$ (counts)")
    plt.xlabel("Horizontal distance $\log_{10}$ (km)")
    plt.ylabel("Temporal distance $\log_{10}$ (s)")
    plt.title("Distribution of horizontal and temporal lags, h=%1.0f km"%(rg))
    plt.savefig("hist_%1.2f.png"%(rg))
    plt.close()
    plt.clf()
        
    ho=h5py.File("hist_%1.2f.h5"%(rg),"w")
    ho["time_diff_log10_s"]=tdist_bins
    ho["horizontal_diff_log10_km"]=hdist_bins
    ho["H"]=H
    ho["hist_hd_x"]=xe
    ho["hist_hd_y"]=ye    
    ho.close()

h.close()
