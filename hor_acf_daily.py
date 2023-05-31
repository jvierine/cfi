import time
import datetime
import numpy as n
import h5py
import matplotlib.pyplot as plt

# our internal modules
import mmaria_read as mr
import cfi_dac as cfi
import cfi_config as c
import mean_wind_est as mw
import os
import stuffr

from mpi4py import MPI
comm=MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()

def avg_hor_acfs(md, # data
                 dcos_thresh=0.8,
                 mean_wind_time_avg=4*3600.0,
                 mean_wind_time=[],
                 h0=90.0,
                 dh=5,
                 ds_h=25.0,
                 dtau=300.0,
                 s_h=n.arange(25.0,500.0,25.0),
                 ds_z=1,
                 t0=0,    # unix seconds for start
                 t1=0,    # unix seconds for end
                 remove_mean=True,
                 name="daily"):


    # what should we do.
    #
    # daily averages of horizontal correlation functions, save them for post processing
    #
    os.system("mkdir -p mpi/%s"%(name))

    n_lags=len(s_h)
    all_acfs=[]
    all_errs=[]
    
    n_avg=0.0
    
    # we need a bit of extra meteors before and after for the mean wind.
    dw=md.read_data(t0=t0-2*3600,t1=t1+2*3600)
    d=md.read_data(t0=t0,t1=t1)    
    n_meas=len(d["t"])
    print("n_meteors %d"%(n_meas))
    if n_meas > 100:
        times,times_h,v,ve,rgs,lat0,lon0,dt,dh,resid=mw.mean_wind_grad(meas=dw,
                                                                       times=mean_wind_time,
                                                                       dt=mean_wind_time_avg,
                                                                       dh=1.0,
                                                                       dz=3.0,
                                                                       max_alt=105,
                                                                       min_alt=78,
                                                                       dcos_thresh=dcos_thresh,
                                                                       ofname="res/tmp-%03d.h5"%(rank),
                                                                       outlier_sigma=4,
                                                                       gradients=False,
                                                                       debug=False,
                                                                       data='dict')
        
            
        meas=cfi.get_meas(meas_file=d,
                          mean_rem=remove_mean,
                          plot_dops=False,
                          dcos_thresh=dcos_thresh,
                          mean_wind_file="res/tmp-%03d.h5"%(rank),
                          data='mmaria')
        
        ih0,idtau,dis_h,acfs,errs,ishs,si_h,names=cfi.hor_acfs(meas,
                                                               h0=h0,
                                                               dh=dh,
                                                               ds_z=ds_z,
                                                               ds_h=ds_h,
                                                               s_h=s_h,
                                                               dtau=dtau,
                                                               LTz=True,
                                                               title=name)

        out_fname="mpi/%s/hacf_res-%06d.h5"%(name,t0)
        os.system("rm -f %s"%(out_fname))        
        
        print("writing %s"%(out_fname))
        ho=h5py.File(out_fname,"w")
        ho["acfs"]=acfs
        ho["t0"]=t0
        ho["t1"]=t1
        ho["errs"]=errs
        ho["s_h"]=si_h
        ho["ds_h"]=ds_h
        ho["ds_z"]=ds_z
        ho["dtau"]=dtau
        ho.close()
    else:
        print("not enough meteors")
            
            
md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()
d0=stuffr.unix2date(b[0])


n_days = int(n.round((b[1]-b[0])/(24*3600)))
t0_floor=stuffr.date2unix(d0.year,d0.month,d0.day,0,0,0)
print("number of days %d floored t0 %1.2f"%(n_days,t0_floor))

# 600 s dtau, 25 km ds_h, 400 km max horizontal lag, 3 day sliding window
name="daily_600s_5day_90km_LTz"

for day_no in range(rank,n_days,size):
    # three day running window with one day increment
    t0 = t0_floor + day_no*24*3600
    t1 = t0_floor + day_no*24*3600 + 5*24*3600

    # 900 second time step for mean wind, but the mean wind averaging window is four hours
    mean_wind_time = n.linspace(t0,t1,num=96*5)
    avg_hor_acfs(md,
                 dcos_thresh=0.8,
                 mean_wind_time_avg=4*3600.0,
                 mean_wind_time=mean_wind_time,
                 h0=90.0,
                 dh=6,
                 ds_h=25.0,            # horizontal distance resolution
                 dtau=600.0,           # lag time resolution
                 s_h=n.arange(0.0,400.0,12.5),
                 ds_z=1,               # vertical lag resolution
                 t0=t0,
                 t1=t1,
                 name=name)


