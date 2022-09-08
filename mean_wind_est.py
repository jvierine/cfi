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
import stuffr

import geoid_const as gc

from mpi4py import MPI
comm=MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()

latdeg2km=gc.latdeg2km #111.321
londeg2km=gc.londeg2km# n.pi*6371.0*n.cos(n.pi*69.0/180.0)/180.0#65.122785 

def mean_wind_grad(meas,
                   times,
                   dt=60*60,
                   t_step=900,
                   dh=1.0,
                   max_alt=105,
                   min_alt=80,
                   dcos_thresh=0.8,
                   gradients=False,
                   ofname="res/mean_wind.h5",
                   min_number_of_measurements=32,
                   debug=False,
                   data='h5file'):  
    
    h=meas
    heights=n.copy(h["heights"])
    heights=heights/1000
    ts=n.copy(h["t"])
    dops=n.copy(h["dops"])
    braggs=n.copy(h["braggs"])
    lats=n.copy(h["lats"])
    lons=n.copy(h["lons"])
    dcoss=n.copy(h["dcos"])
    
    # 1d direction cosine
    dcos2=n.sqrt(dcoss[:,0]**2.0+dcoss[:,1]**2.0)

    # filter by direction cosine
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

#    times=n.arange(int(n.min(ts)),int(n.max(ts)),t_step)
    n_rgs=len(rgs)

    # if no gradients
    n_par=2
    v=n.zeros([n_par,len(times),len(rgs)])
    ve=n.zeros([n_par,len(times),len(rgs)])
    dop_resid=n.zeros([len(times),len(rgs)])
    
    # six unknowns: zonal and meridional mean wind, and their zonal and meridional gradients
    if gradients:
        n_par=6
        v=n.zeros([n_par,len(times),len(rgs)])
        ve=n.zeros([n_par,len(times),len(rgs)])

    # initialize with nan
    ve[:,:,:]=n.nan
    v[:,:,:]=n.nan
    dop_resid[:,:]=n.nan    

    # these measurements are good
    good_meas_idx = n.array([],dtype=n.int)
    bad_meas_idx = n.array([],dtype=n.int)
    
    # for each time step
    for ti,t in enumerate(times):
        print("%d/%d"%(ti,len(times)))
        
        # ridxs is a list of index lists for each range gate        
        ridxs=[]

        # find all measurements that are at the same height and time bin
        for r in rgs:
            hidx=n.where((heights > r)&(heights < (r+dh))&(ts> (t-dt/2))&(ts<(t+dt/2)))[0]
            ridxs.append(hidx)

        # for each height interval
        for ri in range(len(rgs)):
            # ridxs[ri] are indices for measurements that are at this height
            n_meas=len(ridxs[ri])
            if n_meas > min_number_of_measurements:
                A=n.zeros([n_meas,n_par])
            
                # Doppler measurements in units of rad/s
                # 2*pi*f  
                m=-2.0*n.pi*dops[ridxs[ri]]
                # theory
                #
                # v_dop = k \dot (v_0 + (lat-lat0)*v')
                #
                # mean winds
                # Bragg vectors in units of (rad/m)
                # v_u
                A[:,0]=braggs[ridxs[ri],0]
                # v_v
                A[:,1]=braggs[ridxs[ri],1]

                if gradients:
                    #
                    # gradients (how much does velocity change per kilometer
                    # in the zonal and meridional directions
                    #
                    # zon lat grad
                    A[:,2]=braggs[ridxs[ri],0]*(lats[ridxs[ri]]-lat0)*gc.latdeg2km
                    # mer lat grad
                    A[:,3]=braggs[ridxs[ri],1]*(lats[ridxs[ri]]-lat0)*gc.latdeg2km
                    # zon lon grad
                    A[:,4]=braggs[ridxs[ri],0]*(lons[ridxs[ri]]-lon0)*gc.londeg2km
                    # mer lon grad
                    A[:,5]=braggs[ridxs[ri],1]*(lons[ridxs[ri]]-lon0)*gc.londeg2km

                borked=False
                try:
                    xhat=n.linalg.lstsq(A,m)[0]
                except:
                    borked=True

                resid=m-n.dot(A,xhat)
                resid_std=n.median(n.abs(resid))

                # good measurements
                gidx=n.where(n.abs(resid) < 4.0*resid_std)[0]
                bidx=n.where(n.abs(resid) >= 4.0*resid_std)[0]

                good_meas_idx=n.concatenate((good_meas_idx,ridxs[ri][gidx]))
                bad_meas_idx=n.concatenate((bad_meas_idx,ridxs[ri][bidx]))

                # one more iteration with outliers removed
                if len(gidx) > min_number_of_measurements:
                    A2=A[gidx,:]
                    m2=m[gidx]

                    try:
                        xhat2=n.linalg.lstsq(A2,m2)[0]
                    except:
                        borked=True


                    # estimate stdev
                    try:
                        stdev=n.sqrt(n.diag(n.linalg.inv(n.dot(n.transpose(A2),A2))))*resid_std
                    except:
                        borked=True

                    if not borked:
                        resid=m-n.dot(A,xhat2)
                        dop_resid[ti,ri]=n.var(resid)
                        for pi in range(n_par):
                            v[pi,ti,ri]=xhat2[pi]
                            ve[pi,ti,ri]=stdev[pi]

    
        
    times_h=(times-times[0])/3600.0
    dt2=times[1]-times[0]
    dh2=rgs[1]-rgs[0]

    if debug:
        plt.subplot(311)
        plt.pcolormesh(times_h,rgs,n.sqrt(n.transpose(dop_resid)))
        plt.colorbar()
        plt.subplot(312)    
        plt.pcolormesh(times_h,rgs,n.transpose(v[0,:,:]))
        plt.colorbar()    
        plt.subplot(313)        
        plt.pcolormesh(times_h,rgs,n.transpose(v[1,:,:]))
        plt.colorbar()    
        plt.tight_layout()
        plt.show()
    
    ho=h5py.File(ofname,"w")
    ho["times"]=times
    ho["rgs"]=rgs
    ho["v"]=v
    ho["ve"]=ve
    ho["dop_residual"]=dop_resid
    ho["dt"]=dt2
    ho["dh"]=dh2
    ho["lat0"]=lat0
    ho["lon0"]=lon0
    ho["latdeg2km"]=gc.latdeg2km
    ho["londeg2km"]=gc.londeg2km
    ho["good_meas_idx"]=good_meas_idx
    ho["bad_meas_idx"]=bad_meas_idx
    # fraction of bad measurements
    ho["quality"]=float(len(good_meas_idx))/float(len(good_meas_idx)+len(bad_meas_idx))
    ho.close()

    # 
    # This should return a function giving wind as a function of height
    # and time
    #
    return(times,times_h,v,ve,rgs,lat0,lon0,dt2,dh2,dop_resid)

def plot_wind(md,
              t0,
              t1,
              odir="/mnt/data/juha/mmaria_norway/mean_wind",
              dt=900.0,
              gradients=False,
              avg_time=3600*4):
    
    d=md.read_data(t0=t0-avg_time,t1=t1+avg_time)

    n_times=int(n.round((t1-t0)/dt))
    print(n_times)
    times=n.arange(n_times)*dt+t0

    dcos_thresh=0.8
    times,times_h,v,ve,rgs,lat0,lon0,dt,dh,dop_resid=mean_wind_grad(meas=d,
                                                                    times=times,
                                                                    dt=avg_time,
                                                                    dh=1.0,
                                                                    max_alt=105,
                                                                    min_alt=78,
                                                                    ofname="%s/mean_wind-%d.h5"%(odir,t0),
                                                                    gradients=False,
                                                                    dcos_thresh=dcos_thresh)
    v_max=100.0
    dv_max=0.5

    if gradients:    
        plt.figure(figsize=(8,8))
        plt.subplot(321)
        plt.pcolormesh(times_h,rgs,n.transpose(v[0,:,:]),vmin=-v_max,vmax=v_max,cmap="jet")
        plt.title("Zonal mean wind $\overline{u}$ (m/s)")
        plt.xlabel("Time (h)")
        plt.ylabel("Height")
        plt.colorbar()
        plt.subplot(322)
        plt.pcolormesh(times_h,rgs,n.transpose(v[1,:,:]),vmin=-v_max,vmax=v_max,cmap="jet")
        plt.title("Meridional mean wind $\overline{v}$ (m/s)")
        plt.xlabel("Time (h)")
        plt.ylabel("Height")
        plt.colorbar()
        

        plt.subplot(323)
        plt.pcolormesh(times_h,rgs,n.transpose(v[2,:,:]),vmin=-dv_max,vmax=dv_max,cmap="jet")
        plt.title("$d\overline{u}/dy$ (m/s/km)")
        plt.xlabel("Time (h)")
        plt.ylabel("Height")
        plt.colorbar()
        plt.subplot(324)
        
        plt.pcolormesh(times_h,rgs,n.transpose(v[3,:,:]),vmin=-dv_max,vmax=dv_max,cmap="jet")
        plt.title("$d\overline{v}/dy$ (m/s/km)")
        plt.xlabel("Time (h)")
        plt.ylabel("Height")
        plt.colorbar()
        plt.subplot(325)
        plt.pcolormesh(times_h,rgs,n.transpose(v[4,:,:]),vmin=-dv_max,vmax=dv_max,cmap="jet")
        plt.title("$d\overline{u}/dx$ (m/s/km)")
        plt.xlabel("Time (h)")
        plt.ylabel("Height")
        
        plt.colorbar()
        plt.subplot(326)
        plt.pcolormesh(times_h,rgs,n.transpose(v[5,:,:]),vmin=-dv_max,vmax=dv_max,cmap="jet")
        plt.title("$d\overline{v}/dx$ (m/s/km)")
        plt.xlabel("Time (h)")
        plt.ylabel("Height")
        plt.colorbar()        
    else:
        plt.figure(figsize=(8,8))
        plt.subplot(311)
        plt.pcolormesh(times_h,rgs,n.transpose(v[0,:,:]),vmin=-v_max,vmax=v_max,cmap="jet")
        plt.title("Zonal mean wind $\overline{u}$ (m/s)")
        plt.xlabel("Time (h since %s)"%(stuffr.unix2datestr(t0)))
        plt.ylabel("Height")
        plt.colorbar()
        
        plt.subplot(312)
        plt.pcolormesh(times_h,rgs,n.transpose(v[1,:,:]),vmin=-v_max,vmax=v_max,cmap="jet")
        plt.title("Meridional mean wind $\overline{v}$ (m/s)")
        plt.xlabel("Time (h)")
        plt.ylabel("Height")
        plt.colorbar()

        plt.subplot(313)
        plt.pcolormesh(times_h,rgs,n.transpose(n.sqrt(dop_resid)),vmin=0,vmax=100.0,cmap="plasma")
        plt.title("Doppler residual RMS (rad/s)")
        plt.xlabel("Time (h)")
        plt.ylabel("Height")
        plt.colorbar()
        

    plt.tight_layout()
    plt.savefig("%s/mean_%d.png"%(odir,t0))
    plt.clf()
    plt.close()

def plot_all():
    md=mr.mmaria_data(c.data_directory)#for many files in a directory
    b=md.get_bounds()
    
    n_days=int(n.round((b[1]-b[0])/(24*3600)))
    
    bt0=int(n.round(b[0]/(24*3600)))*24*3600
    
    for di in range(rank,n_days,size):
        plot_wind(md,
                  t0=bt0+24*3600*di,
                  t1=bt0+24*3600*(di+1),
                  odir="/mnt/data/juha/mmaria_norway/mean_wind",
                  avg_time=3600*4)

if __name__ == "__main__":
    plot_all()
