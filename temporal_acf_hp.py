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
                      dstep=4,
                      mean_rem=True,
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
#            dstep=n.max([1,int(n.floor(0.5*max_lag/(24*3600.0)))])
            for day in n.arange(1,31,dstep):
                pars.append([year,month,int(day)])

            
    n_pars=len(pars)
    print("n_runs %d"%(n_pars))    
    
    for pi in range(rank,n_pars,size):
        year=pars[pi][0]
        month=pars[pi][1]
        day=pars[pi][2]
        d0=datetime.date(year,month,day)
        t0=time.mktime(d0.timetuple())
    
        d=md.read_data(t0=t0,t1=t0+max_lag+dstep*24*3600)
        n_meas=len(d["t"])
        print("n_meteors %d"%(n_meas))
        if n_meas > 100:

            if mean_rem:
                mw_fname="res/mw-%1.6f.h5"%(n.random.rand(1))
                n_mw_t=int(n.ceil((max_lag+dstep*24*3600)/900.0))
                mw.mean_wind_grad(d,
                                  times=n.linspace(t0,t0+max_lag+dstep*24*3600,num=n_mw_t),
                                  dt=4*3600.0,
                                  t_step=900,
                                  dh=1.0,
                                  max_alt=105,
                                  min_alt=80,
                                  dcos_thresh=0.8,
                                  ofname=mw_fname,
                                  min_number_of_measurements=32)

            meas=cfi.get_meas(meas_file=d,
                              mean_rem=mean_rem,
                              plot_dops=False,
                              dcos_thresh=dcos_thresh,
                              mean_wind_file=mw_fname,
                              data='mmaria')
            if mean_rem:
                os.system("rm %s"%(mw_fname))
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

def summer_hp():
    avg_temporal_acfs(dcos_thresh=0.7,
                      h0=90.0,
                      dh=2,
                      ds_h=50.0,
                      dtau=60.0,
                      tau=n.arange(0.0,8*3600,60.0),
                      ds_z=1,
                      years=[2018,2019],
                      months=[5,6,7],
                      dstep=7,
                      mean_rem=True,
                      name="summer_tacf_8_60_50km")
    

#winter_small_scale()
summer_hp()
