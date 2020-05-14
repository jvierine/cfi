#
# 2d correlation function
#
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

def avg_mer_acfs(md, # data
                 dcos_thresh=0.8,
                 h0=90.0,
                 dh=2,
                 ds_x=100.0,
                 ds_y=100.0,
                 dtau=600.0,
                 s_y=n.arange(-500,500.0,100.0),
                 tau=n.arange(-3600,3600,600),
                 ds_z=1,
                 years=[2018,2019,2020],                 
                 months=[5,6,7],
                 name="summer_2d_mer_time"):

    if rank == 0:
        os.system("rm mpi/%s/*.h5"%(name))        
        os.system("mkdir -p mpi/%s"%(name))
    comm.Barrier()

    n_mer_lags=len(s_y)
    n_time_lags=len(tau)
    
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
        
        d=md.read_data(t0=t0,t1=t0+3*31*24*3600.0)
        n_meas=len(d["t"])
        print("n_meteors %d"%(n_meas))
        
        if n_meas > 100:
            meas=cfi.get_meas(meas_file=d,
                              mean_rem=False,
                              plot_dops=False,
                              dcos_thresh=dcos_thresh,
                              mean_wind_file="res/tmp-%03d.h5"%(rank),
                              data='mmaria')
            
            acfs=n.zeros([6,n_mer_lags,n_time_lags])
            errs=n.zeros([6,n_mer_lags,n_time_lags])
            
            for ti,tlag in enumerate(tau):
                for yi,y_lag in enumerate(s_y):
                    print("t %1.2f mer %1.2f"%(tlag,y_lag))
                    acf, err, itaus, ixs, iys, izs, is_h=cfi.cfi(meas,
                                                                 h0=h0,
                                                                 dh=2,
                                                                 ds_z=1.0,
                                                                 ds_x=ds_x,
                                                                 ds_y=ds_y,
                                                                 dtau=dtau,
                                                                 tau=tlag,
                                                                 s_x=0.0,
                                                                 s_y=y_lag,
                                                                 horizontal_dist=False)

                    print(acf)
                    acfs[:,yi,ti]=acf
                    errs[:,yi,ti]=err
                    
            ho=h5py.File("mpi/%s/2d_mer_time_res-%06d.h5"%(name,pi),"w")
            ho["acfs"]=acfs
            ho["errs"]=errs
            ho["s_y"]=s_y
            ho["tau"]=tau
            ho["ds_x"]=ds_x
            ho["ds_y"]=ds_y
            ho["dtau"]=dtau
            ho.close()
            
md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()


def summer():
    avg_mer_acfs(md, 
                 dcos_thresh=0.8,
                 h0=90.0,
                 dh=2,
                 ds_x=100.0,
                 ds_y=100.0,
                 dtau=900.0,
                 s_y=n.arange(-400.0, 400.0, 50.0),
                 tau=n.arange(-12*3600,12*3600,900.0),
                 ds_z=1,
                 years=[2018,2019,2020],
                 months=[5],
                 name="summer_2d_mer_time_long")

summer()
