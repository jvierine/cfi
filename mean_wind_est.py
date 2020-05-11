#!/usr/bin/env python
#
# simple mean horizontal wind
#
import h5py
import matplotlib.pyplot as plt
import glob
from mpl_toolkits.mplot3d import Axes3D
import numpy as n

import mmaria_read as mr
import cfi_dac as cfi
import cfi_config as c
import time
import datetime

import geoid_const as gc

latdeg2km=gc.latdeg2km #111.321
londeg2km=gc.londeg2km# n.pi*6371.0*n.cos(n.pi*69.0/180.0)/180.0#65.122785 

# high time resolution mean wind
#t_avg=1800
# low time resolution mean wind (used for high pass filtering)
#lf_t_avg=4*3600.0
#dcos_thresh=0.8
#ofname="res/mean_wind_4h.h5"

def mean_wind(meas="res/simone_nov2018_multilink_juha_30min_1000m.h5",
              dt=60*60,t_step=900,dh=1.0,max_alt=105,min_alt=80,dcos_thresh=0.8,
              ofname="res/mean_wind.h5",
              data='h5file'):   #data is the type of input values.
    
    if data=='h5file':
        h=h5py.File(meas,"r")
        print(h.keys())
        heights=n.copy(h["heights"].value)
        ts=n.copy(h["t"].value)
        dops=n.copy(h["dops"].value)
        braggs=n.copy(h["braggs"].value)
        lats=n.copy(h["lats"].value)
        if "dcos" in h.keys():
            dcoss=n.copy(h["dcos"].value)
        else:
            dcoss=n.zeros([len(heights),2])
        
        lons=n.copy(h["lons"].value)
    
    else:   #for now only for mmaria_read.py
        h=meas
        heights=n.copy(h["heights"])
        heights=heights/1000
        ts=n.copy(h["t"])
        dops=n.copy(h["dops"])
        braggs=n.copy(h["braggs"])
        lats=n.copy(h["lats"])
        if "dcos" in h.keys():
            dcoss=n.copy(h["dcos"])
        else:
            dcoss=n.zeros([len(heights),2])
        
        lons=n.copy(h["lons"])
        
        
        
    dcos2=n.sqrt(dcoss[:,0]**2.0+dcoss[:,1]**2.0)
    
    ok_idx=n.where(dcos2 < dcos_thresh)[0]
    ts=ts[ok_idx]
    lats=lats[ok_idx]
    lons=lons[ok_idx]
    heights=heights[ok_idx]
    braggs=braggs[ok_idx,:]
    dops=dops[ok_idx]
    dcoss=dcoss[ok_idx,:]

    lat0=n.median(lats)
    lon0=n.median(lons)
    
    rgs=n.arange(min_alt,max_alt,int(dh))
    times=n.arange(int(n.min(ts)),int(n.max(ts)),t_step)
    n_rgs=len(rgs)
    
    np=2
    
    v=n.zeros([np,len(times),len(rgs)])
    ve=n.zeros([np,len(times),len(rgs)])
    ve[:,:,:]=n.nan
    v[:,:,:]=n.nan

    for ti,t in enumerate(times):
        print("%d/%d"%(ti,len(times)))
        ridxs=[]
        n_r=[]
        for r in rgs:
            hidx=n.where((heights > r)&(heights < (r+dh))&(ts>t-dt/2)&(ts<(t+dt/2)))[0]
            ridxs.append(hidx)
            n_r.append(len(hidx))
            
        n_r=n.array(n_r)
        n_meas=n.sum(n_r)
        # v_x = v_x
        # v_y = v_y
        A=n.zeros([n_meas,np*n_rgs])
        m=n.zeros(n_meas)
        n_m = 0
        g_ridx2=[]
        g_ridx=[]
        
        for ri in range(len(rgs)):
            m_idx=n.arange(n_m,n_m+len(ridxs[ri]))
            # meas
            m[m_idx]=-2.0*n.pi*dops[ridxs[ri]]
            # theory
            # v_u
            A[m_idx,ri*np+0]=braggs[ridxs[ri],0]
            # v_v
            A[m_idx,ri*np+1]=braggs[ridxs[ri],1]
        
            n_m+=len(ridxs[ri])
            if len(ridxs[ri]) > 10:
                for pi in range(np):
                    g_ridx2.append(ri*np+pi)
                g_ridx.append(ri)
            else:
                pass
#                print("bad range %d"%(ri))
        g_ridx2=n.array(g_ridx2)
        g_ridx=n.array(g_ridx)
    
        if len(g_ridx2) > 0:
            xhat=n.linalg.lstsq(A[:,g_ridx2],m)[0]
            resid=m-n.dot(A[:,g_ridx2],xhat)
            gidx=n.where(n.abs(resid) < 100.0)[0]
            #        plt.plot(resid)
            #       plt.show()
            
            g_ridx_f=[]
            g_ridx2_f=[]        
            for ri in g_ridx:
                n_meas_per_rg=len(n.where(n.abs(A[:,np*ri])>0)[0])
                if n_meas_per_rg > np:
                    g_ridx_f.append(ri)                
                for pi in range(np):
                    g_ridx2_f.append(np*ri+pi)

            g_ridx_f=n.array(g_ridx_f,dtype=n.int)
            g_ridx2_f=n.array(g_ridx2_f,dtype=n.int)

            # estimate stdev
            stdev=n.sqrt(n.diag(n.linalg.inv(n.dot(n.transpose(A[:,g_ridx2_f]),A[:,g_ridx2_f]))))*n.std(resid)

            n_gridx=len(g_ridx_f)
            for pi in range(np):
                v[pi,ti,g_ridx_f]=xhat[np*n.arange(n_gridx)+pi]
                ve[pi,ti,g_ridx_f]=stdev[np*n.arange(n_gridx)+pi]

    times_h=(times-times[0])/3600.0
    dt2=times[1]-times[0]
    dh2=rgs[1]-rgs[0]

    ho=h5py.File(ofname,"w")
    ho["times"]=times
    ho["rgs"]=rgs
    ho["v"]=v
    #ho["v_fluct"]=v2-v
    ho["ve"]=ve
    ho["dt"]=dt2
    ho["dh"]=dh2
    ho["lat0"]=lat0
    ho["lon0"]=lon0
    ho.close()

    
    return(times,times_h,v,ve,rgs,lat0,lon0,dt2,dh2)


if __name__ == "__main__":

    md=mr.mmaria_data(c.data_directory)#for many files in a directory
    b=md.get_bounds()

    d=md.read_data_date(d0=datetime.date(2019,1,5),d1=datetime.date(2019,1,9))

    lf_t_avg=4*3600
    t_avg=900.0
    times,times_h,v,ve,rgs,lat0,lon0,dt,dh=mean_wind(meas=d, dt=lf_t_avg,dh=1.0,max_alt=105,min_alt=78,dcos_thresh=dcos_thresh,data='dict')
    #times,times_h,v,ve,rgs,lat0,lon0,dt,dh=mean_wind(dt=1800,dh=1.0,max_alt=105,min_alt=78,dcos_thresh=0.8)
    times2,times_h2,v2,ve2,rgs2,lat02,lon02,dt2,dh2=mean_wind(meas=d, dt=t_avg,dh=1.0,max_alt=105,min_alt=78,dcos_thresh=dcos_thresh,data='dict')

    print(dh)


    plt.figure(figsize=(8,8))

    plt.subplot(322)
    plt.pcolormesh(times_h,rgs,n.transpose(v2[0,:,:]),vmin=-70,vmax=70,cmap="jet")
    plt.title("Zonal mean wind (m/s)")
    #plt.xlabel("Time (h)")
    plt.ylabel("%d minute average"%(t_avg/60.0))
    plt.colorbar()

    plt.subplot(321)
    plt.pcolormesh(times_h,rgs,n.transpose(v2[1,:,:]),vmin=-70,vmax=70,cmap="jet")
    plt.title("Meridional mean wind (m/s)")
    #plt.xlabel("Time (h)")
    plt.ylabel("Altitude (km)")
    plt.colorbar()

    plt.subplot(324)
    plt.pcolormesh(times_h,rgs,n.transpose(v[0,:,:]),vmin=-70,vmax=70,cmap="jet")
    #plt.title("Zonal mean wind (m/s)")
    #plt.xlabel("Time (h)")
    plt.ylabel("4 hour average")
    plt.colorbar()

    plt.subplot(323)
    plt.pcolormesh(times_h,rgs,n.transpose(v[1,:,:]),vmin=-70,vmax=70,cmap="jet")
    #plt.title("Meridional mean wind (m/s)")
    #plt.xlabel("Time (h)")
    plt.ylabel("Altitude (km)")
    plt.colorbar()

    plt.subplot(326)
    plt.pcolormesh(times_h,rgs,n.transpose(v2[0,:,:]-v[0,:,:]),vmin=-25,vmax=25,cmap="jet")
    #plt.title("Zonal mean wind (m/s)")
    plt.xlabel("Time (h)")
    plt.ylabel("Residual")
    plt.colorbar()

    plt.subplot(325)
    plt.pcolormesh(times_h,rgs,n.transpose(v2[1,:,:]-v[1,:,:]),vmin=-25,vmax=25,cmap="jet")
    #plt.title("Meridional mean wind (m/s)")
    plt.xlabel("Time (h)")
    #plt.ylabel("Residual")
    plt.ylabel("Altitude (km)")
    plt.colorbar()

    plt.tight_layout()
    plt.savefig("mean_wind.png")
    plt.show()
    

