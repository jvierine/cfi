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
    os.system("mkdir -p %s/%s/%s"%(c.data_directory,c.data_prefix,name))
    os.system("mkdir -p %s/%s/%s/tmp"%(c.data_directory,c.data_prefix,name))    

    n_lags=len(s_h)
    all_acfs=[]
    all_errs=[]
    
    n_avg=0.0

    dw=None
    if remove_mean:
        # we need a bit of extra meteors before and after for the mean wind.
        dw=md.read_data(t0=t0-2*3600,t1=t1+2*3600)

        if False:
            # a little debug plot to show what meteors are read
            plt.subplot(121)
            plt.plot((dw["t"]-t0)/24/3600.0,dw["heights"],".")
            plt.subplot(122)
            plt.plot(dw["lats"],dw["lons"],".")
            plt.tight_layout()
            plt.show()
    print(dw)
    d=md.read_data(t0=t0,t1=t1)    
    n_meas=len(d["t"])
    print("n_meteors %d"%(n_meas))
    mw_fname="%s/%s/%s/tmp/tmp-%05d.h5"%(c.data_directory,c.data_prefix,name,rank)
    if n_meas > 100:
        if remove_mean:
            mwr=mw.mean_wind_grad(meas=dw,
                                  times=mean_wind_time,
                                  dt=mean_wind_time_avg,
                                  dh=1.0,
                                  dz=3.0,
                                  max_alt=110,
                                  min_alt=70,
                                  dcos_thresh=dcos_thresh,
                                  ofname=mw_fname,
                                  outlier_sigma=4,    # don't be overly strict about removing measurements that don't conform with the mean wind
                                  gradients=False,
                                  debug=False)

        #  note that outlier removal is done in cfi.get_meas
        meas=cfi.get_meas(meas_file=d,
                          mean_rem=remove_mean,
                          plot_dops=False,
                          dcos_thresh=dcos_thresh,
                          mean_wind_file=mw_fname,
                          data='mmaria')
        
        ih0,idtau,dis_h,acfs,errs,ishs,si_h,names=cfi.hor_acfs(meas,
                                                               h0=h0,
                                                               dh=dh,
                                                               ds_z=ds_z,
                                                               ds_h=ds_h,
                                                               s_h=s_h,
                                                               dtau=dtau,
                                                               LTz=c.horizontal_correlation_LTz,
                                                               title=name)

        out_fname="%s/%s/%s/hacf_res-%06d.h5"%(c.data_directory,c.data_prefix,name,t0)
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
        try:
            os.system("rm %s"%(mw_fname))
        except:
            print("couldn't remove %s"%(mw_fname))
            pass
    else:
        print("not enough meteors")
            
            
md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()
d0=stuffr.unix2date(b[0])

n_days = int(n.round((b[1]-b[0])/(24*3600)))
t0_floor=stuffr.date2unix(d0.year,d0.month,d0.day,0,0,0)
print("number of days %d floored t0 %1.2f"%(n_days,t0_floor))

for h0 in c.horizontal_correlation_heights:
    # 600 s dtau, 25 km ds_h, 400 km max horizontal lag, 3 day sliding window
    #name="daily_600s_7day_82km_LTz"
    name="hcor_%1.2f"%(h0)
        
    for day_no in range(rank,n_days,size):
        # three day running window with one day increment
        t0 = t0_floor + day_no*24*3600
        t1 = t0_floor + day_no*24*3600 + c.horizontal_correlation_n_days*24*3600
        
        # 900 second time step for mean wind, but the mean wind averaging window is four hours
        mean_wind_time = n.linspace(t0,t1,num=int(96*c.horizontal_correlation_n_days))
        print("averaging correlation function of %d days"%(c.horizontal_correlation_n_days))
        try:
            avg_hor_acfs(md,
                         dcos_thresh=c.dcos_thresh,
                         mean_wind_time_avg=c.horizontal_correlation_mean_wind_avg_time*3600.0,
                         mean_wind_time=mean_wind_time,
                         h0=h0,
                         dh=c.horizontal_correlation_dh, # look for meteors between [h0-dh/2,h0+dh/2]
                         ds_h=25.0,            # horizontal distance resolution
                         dtau=c.horizontal_correlation_dtau,           # lag time resolution
                         s_h=n.arange(0.0,400.0,12.5),
                         ds_z=1,               # vertical lag resolution
                         t0=t0,
                         t1=t1,
                         remove_mean=c.high_pass_filter,
                         name=name)
        except:
            print("couldn't calculate correlation function!")
            pass



