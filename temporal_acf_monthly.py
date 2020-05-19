import time
import datetime
import numpy as n
import h5py
import matplotlib.pyplot as plt
import os

# our internal modules
import mmaria_read as mr
import cfi_dac as cfi
import cfi_config as c
import mean_wind_est as mw

from mpi4py import MPI
comm=MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()

md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()



def avg_temporal_acfs(dcos_thresh=0.8,
                      mean_wind_time_avg=4*3600.0,
                      h0=90.0,
                      dh=5,
                      ds_h=50.0,
                      dtau=300.0,
                      tau=n.arange(0.0,7*24*3600,300.0),
                      ds_z=1,
                      years=[2018,2019,2020],
                      months=[5,6,7],
                      name="summer_tacf"):
    if rank == 0:
#        os.system("rm mpi/%s/tacf_res*.h5"%(name))
        os.system("mkdir -p mpi/%s"%(name))
    comm.Barrier()

    n_lags=len(tau)
    
    max_lag=n.max(tau)

    pars=[]
    for year in years:
        for month in months:
            pars.append([year,month,1])

            
    n_pars=len(pars)
    print("n_runs %d"%(n_pars))    
    
    for pi in range(rank,n_pars,size):
        year=pars[pi][0]
        month=pars[pi][1]
        day=pars[pi][2]
        d0=datetime.date(year,month,day)
        t0=time.mktime(d0.timetuple())
    
        d=md.read_data(t0=t0,t1=t0+max_lag+31*24*3600)
        n_meas=len(d["t"])
        print("n_meteors %d"%(n_meas))
        if n_meas > 100:
            meas=cfi.get_meas(meas_file=d,
                              mean_rem=False,
                              plot_dops=False,
                              dcos_thresh=dcos_thresh,
                              mean_wind_file="res/tmp.h5",
                              data='mmaria')
                    
            acfs,errs,tau,dtau,ds_h,names=cfi.temporal_acfs(meas,
                                                            h0=h0,
                                                            tau=tau,
                                                            dtau=dtau,
                                                            min_dt=10,
                                                            ds_h=ds_h)
            ofname="mpi/%s/tacf_res-%06d.h5"%(name,pi)
            print(ofname)
            ho=h5py.File(ofname,"w")
            ho["acfs"]=acfs
            ho["errs"]=errs
            ho["tau"]=tau
            ho["dtau"]=dtau
            ho["ds_h"]=ds_h
            ho["ds_z"]=ds_z
            ho["h0"]=h0
            ho["dh"]=dh
            ho.close()

def summer_day():
    avg_temporal_acfs(dcos_thresh=0.8,
                      h0=90.0,
                      dh=2,
                      ds_h=50.0,
                      dtau=100.0,
                      tau=n.arange(0.0,24*3600,100.0),
                      ds_z=1,
                      years=[2018,2019],
                      months=[5,6,7],
                      dstep=2,
                      name="summer_tacf_1_120_50km")
    
def summer_small_scale():
    avg_temporal_acfs(dcos_thresh=0.8,
                      h0=90.0,
                      dh=2,
                      ds_h=100.0,
                      dtau=1800.0,
                      tau=n.arange(0.0,14*24*3600,1800.0),
                      ds_z=1,
                      years=[2018,2019],
                      months=[5,6,7],
                      name="summer_tacf_koki")
    
def winter_small_scale():
    avg_temporal_acfs(dcos_thresh=0.8,
                      h0=90.0,
                      dh=2,
                      ds_h=100.0,
                      dtau=1800.0,
                      tau=n.arange(0.0,14*24*3600,1800.0),
                      ds_z=1,
                      years=[2018,2019,2020],
                      months=[11,12,1],
                      name="winter_tacf_koki")
    

#winter_small_scale()
summer_small_scale()
