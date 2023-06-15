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
import os
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
                   dh=1.0,   # vertical sample spacing
                   dz=3.0,   # maximum deviation in vertical direction
                   max_alt=105,
                   min_alt=80,
                   dcos_thresh=0.8,
                   outlier_sigma=4,
                   gradients=False,
                   ofname="res/mean_wind.h5",
                   min_number_of_measurements=10,
                   debug=False):

    
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
#        print("%d/%d"%(ti,len(times)))
        
        # ridxs is a list of index lists for each range gate        
        ridxs=[]

        # find all measurements that are at the same height and time bin
        for r in rgs:
            # the range gate is centered at range r and has a width of dh
            # the time step is centered at t and has a width of dt
            hidx=n.where((heights > (r-dz/2.0))&(heights < (r+dz/2.0))&(ts> (t-dt/2))&(ts<(t+dt/2)))[0]
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
                # robust estimator for standard deviation
                resid_std=0.7*n.median(n.abs(resid))

                # good measurements
                # outlier rejection, when making the estimate of mean wind
                # by default, 100 standard deviations is an outlier
                gidx=n.where(n.abs(resid) < outlier_sigma*resid_std)[0]
                bidx=n.where(n.abs(resid) >= outlier_sigma*resid_std)[0]

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
        plt.title("Dopppler residual stdev (Hz)")
        plt.colorbar()
        plt.subplot(312)    
        plt.pcolormesh(times_h,rgs,n.transpose(v[0,:,:]))
        plt.title("Zonal wind (m/s)")        
        plt.colorbar()    
        plt.subplot(313)        
        plt.pcolormesh(times_h,rgs,n.transpose(v[1,:,:]))
        plt.title("Meridional wind (m/s)")                
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
    ho["good_meas_idx"]=ok_idx[good_meas_idx]
    ho["bad_meas_idx"]=ok_idx[bad_meas_idx]
    # fraction of bad measurements
    ho["quality"]=float(len(good_meas_idx))/float(len(good_meas_idx)+len(bad_meas_idx))
    ho.close()

    # only use good measurements
    good_meas_idx=n.unique(good_meas_idx)
    bad_meas_idx=n.unique(bad_meas_idx)
    # remove measurements that where considered bad at any possible point of doing the mean wind estimates
    good_meas_idx = n.setdiff1d(good_meas_idx,bad_meas_idx)
    # 
    # This should return a function giving wind as a function of height
    # and time
    #
    #    return(times,times_h,v,ve,rgs,lat0,lon0,dt2,dh2,dop_resid)
    return({"times":times,
            "times_h":times_h,
            "v":v,
            "ve":ve,
            "rgs":rgs,
            "lat0":lat0,
            "lon0":lon0,
            "dt2":dt2,
            "dh2":dh2,
            "dop_resid":dop_resid,
            "good_idx":ok_idx[good_meas_idx],
            "bad_idx":ok_idx[bad_meas_idx]})

def plot_wind(md,
              t0,
              t1,
              odir="/mnt/data/juha/mmaria_norway/mean_wind",
              dt=900.0,
              gradients=True,
              avg_time=3600*4):
    
    d=md.read_data(t0=t0-avg_time,t1=t1+avg_time)
    
    n_times=int(n.round((t1-t0)/dt))
    print(n_times)
    times=n.arange(n_times)*dt+t0

    dcos_thresh=0.8
    #times,times_h,v,ve,rgs,lat0,lon0,dt,dh,dop_resid
    mwr=mean_wind_grad(meas=d,
                        times=times,
                        dt=avg_time,
                        dh=1.0,
                        max_alt=105,
                        min_alt=78,
                        ofname="%s/mean_wind-%d.h5"%(odir,t0),
                        gradients=True,
                        dcos_thresh=dcos_thresh)
    times=mwr["times"]
    times_h=mwr["times_h"]
    v=mwr["v"]
    ve=mwr["ve"]
    rgs=mwr["rgs"]
    lat0=mwr["lat0"]
    lon0=mwr["lon0"]
    dt=mwr["dt"]
    dh=mwr["dh"]
    dop_resid=mwr["dop_resid"]
    
    v_max=200.0
    dv_max=0.5


    # zonal
    qu=n.copy(n.transpose(v[0,:,:]))
    # meridional
    qv=n.copy(n.transpose(v[1,:,:]))
    qu[qu<-v_max]=-v_max
    qv[qv<-v_max]=-v_max
    qu[qu>v_max]=v_max
    qv[qv>v_max]=v_max
    
    plt.quiver(times_h, rgs, qu, qv)
    plt.title("Horizontal wind (right is east)")
    plt.xlabel("Time (UT hour)")
    plt.ylabel("Height (km)")
    plt.tight_layout()
    plt.savefig("%s/quiv_%d.png"%(odir,t0))
    plt.clf()
    plt.close()

    
    if gradients:    
        plt.figure(figsize=(20,20))
        plt.subplot(321)
        plt.pcolormesh(times_h,rgs,n.transpose(v[0,:,:]),vmin=-v_max,vmax=v_max,cmap="jet")
        plt.title("Zonal mean wind $\overline{u}$ (m/s)")
        plt.xlabel("Time (h since %s)"%(stuffr.unix2datestr(t0)))        
#        plt.xlabel("Time (h)")
        plt.ylabel("Height (km)")
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
        plt.subplot(221)
        plt.pcolormesh(times_h,rgs,n.transpose(v[0,:,:]),vmin=-v_max,vmax=v_max,cmap="jet")
        plt.title("Zonal mean wind $\overline{u}$ (m/s)")
        plt.xlabel("Time (h since %s)"%(stuffr.unix2datestr(t0)))
        plt.ylabel("Height")
        plt.colorbar()
        
        plt.subplot(222)
        plt.pcolormesh(times_h,rgs,n.transpose(v[1,:,:]),vmin=-v_max,vmax=v_max,cmap="jet")
        plt.title("Meridional mean wind $\overline{v}$ (m/s)")
        plt.xlabel("Time (h)")
        plt.ylabel("Height")
        plt.colorbar()

        plt.subplot(223)
        plt.pcolormesh(times_h,rgs,n.transpose(n.sqrt(dop_resid)),vmin=0,vmax=50.0,cmap="plasma")
        plt.title("Doppler residual RMS (rad/s))")
        plt.xlabel("Time (h)")
        plt.ylabel("Height")
        plt.colorbar()
        
        plt.subplot(224)
        plt.pcolormesh(times_h,rgs,n.sqrt(n.transpose(v[1,:,:])**2.0 + n.transpose(v[0,:,:])**2.0),vmin=0,vmax=v_max,cmap="plasma")
        plt.title("Velocity magnitude (m/s)")
        plt.xlabel("Time (h)")
        plt.ylabel("Height")
        plt.colorbar()
        

    plt.tight_layout()
    plt.savefig("%s/mean_%d.png"%(odir,t0))
    plt.clf()
    plt.close()

    plt.figure(figsize=(8,8))
    plt.subplot(221)
    # v = c*df/2/f_rad
    plt.pcolormesh(times_h,rgs,n.sqrt(n.transpose(dop_resid)),vmin=0,vmax=50,cmap="jet")
    plt.title("Dopppler residual stdev (rad/s)")
    plt.colorbar()
    plt.subplot(222)    
    plt.pcolormesh(times_h,rgs,n.transpose(v[0,:,:]),vmin=-v_max,vmax=v_max,cmap="jet")
    plt.title("Zonal wind (m/s)")        
    plt.colorbar()    
    plt.subplot(223)        
    plt.pcolormesh(times_h,rgs,n.transpose(v[1,:,:]),vmin=-v_max,vmax=v_max,cmap="jet")
    plt.title("Meridional wind (m/s)")                
    plt.colorbar()    


    plt.subplot(224)        
    plt.pcolormesh(times_h,rgs,n.sqrt(n.transpose(dop_resid))/(n.sqrt(n.transpose(v[1,:,:])**2.0+n.transpose(v[0,:,:])**2.0)),vmin=0,vmax=3,cmap="jet")
    plt.title("Normalized Doppler res stdev ((rad/s)/(m/s))")
    plt.colorbar()    
    
    plt.tight_layout()
    plt.savefig("%s/resid_%d.png"%(odir,t0))
    plt.clf()
    plt.close()
    
def plot_all():
    md=mr.mmaria_data(c.data_directory)#for many files in a directory
    b=md.get_bounds()
    
    n_days=int(n.round((b[1]-b[0])/(24*3600)))
    
    bt0=int(n.round(b[0]/(24*3600)))*24*3600

    os.system("mkdir -p %s/mean_wind"%(c.plot_directory))
    for di in range(rank,n_days,size):
        try:
            plot_wind(md,
                      t0=bt0+24*3600*di,
                      t1=bt0+24*3600*(di+1),
                      odir="%s/mean_wind"%(c.plot_directory),
                      dt=600,
                      avg_time=3600)
        except:
            print("error with %s"%(di))

if __name__ == "__main__":
    plot_all()
