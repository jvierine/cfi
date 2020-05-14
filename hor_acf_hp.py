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

from mpi4py import MPI
comm=MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()

def avg_hor_acfs(md, # data
                 dcos_thresh=0.8,
                 mean_wind_time_avg=4*3600.0,
                 h0=90.0,
                 dh=5,
                 ds_h=25.0,
                 dtau=300.0,
                 s_h=n.arange(25.0,500.0,25.0),
                 ds_z=1,
                 years=[2018,2019,2020],                 
                 months=[5,6,7],
                 n_days=7,
                 avg_time=4*3600.0,
                 name="summer_hacf_hp",
                 mean_rem=False):

    if rank == 0:
        os.system("rm mpi/%s/*.h5"%(name))        
        os.system("mkdir -p mpi/%s"%(name))
    comm.Barrier()

    n_lags=len(s_h)
    all_acfs=[]
    all_errs=[]
    
    n_avg=0.0


    pars=[]
    for year in years:
        for month in months:
            for day in range(0,31,n_days):
                pars.append([year,month,day+1])
    n_pars=len(pars)
    print("n_runs %d"%(n_pars))

    for pi in range(rank,n_pars,size):
        year=pars[pi][0]
        month=pars[pi][1]
        day=pars[pi][2]
        d0=datetime.date(year,month,day)
        t0=time.mktime(d0.timetuple())
        
        d=md.read_data(t0=t0,t1=t0+n_days*24*3600.0)
        n_meas=len(d["t"])
        print("n_meteors %d"%(n_meas))
        if n_meas > 100:
            mw_fname="res/mw-%1.6f.h5"%(n.random.rand(1))
            if mean_rem:
                n_mw_t=int(n.ceil(n_days*24*3600/900.0))
                mw.mean_wind_grad(d,
                                  times=n.linspace(t0,t0+n_days*24*3600,num=n_mw_t),
                                  dt=avg_time,
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
            

            ih0,idtau,dis_h,acfs,errs,ishs,si_h,names=cfi.hor_acfs(meas,
                                                                   h0=h0,
                                                                   dh=dh,
                                                                   ds_z=ds_z,
                                                                   ds_h=ds_h,
                                                                   s_h=s_h,
                                                                   dtau=dtau,
                                                                   title=name)
            
            ho=h5py.File("mpi/%s/hacf_res-%06d.h5"%(name,pi),"w")
            ho["acfs"]=acfs
            ho["errs"]=errs
            ho["s_h"]=si_h
            ho["ds_h"]=ds_h
            ho["ds_z"]=ds_z
            ho["dtau"]=dtau
            ho.close()
            
md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()

def summer_hp():
    avg_hor_acfs(md,
                 dcos_thresh=0.8,
                 mean_wind_time_avg=4*3600.0,
                 h0=90.0,
                 dh=2,
                 ds_h=25.0,
                 dtau=900.0,
                 s_h=n.arange(0.0,400.0,12.5)+12.5,
                 ds_z=1,
                 years=[2018,2019],
                 months=[5,6,7],
                 n_days=7,
                 name="summer_hp_hacf",
                 mean_rem=True)

summer_hp()
