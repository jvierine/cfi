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

def avg_mer_zon_acfs(md, # data
                     dcos_thresh=0.8,
                     h0=90.0,
                     dh=2,
                     ds_x=100.0,
                     ds_y=100.0,
                     dtau=600.0,
                     s_x=n.arange(-500,500.0,100.0),
                     s_y=n.arange(-500,500.0,100.0),
                     ds_z=1,
                     years=[2018,2019,2020],                 
                     months=[5,6,7],
                     name="summer_2d_mer_zon"):

    if rank == 0:
        os.system("rm mpi/%s/*.h5"%(name))        
        os.system("mkdir -p mpi/%s"%(name))
    comm.Barrier()

    n_mer_lags=len(s_y)
    n_zon_lags=len(s_x)
    
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
            
            acfs=n.zeros([6,n_mer_lags,n_zon_lags])
            errs=n.zeros([6,n_mer_lags,n_zon_lags])
            
            for yi,y_lag in enumerate(s_y):
                for xi,x_lag in enumerate(s_x):
                    print("zon %1.2f mer %1.2f"%(y_lag,x_lag))
                    acf, err, itaus, ixs, iys, izs, is_h=cfi.cfi(meas,
                                                                 h0=h0,
                                                                 dh=2,
                                                                 ds_z=1.0,
                                                                 ds_x=ds_x,
                                                                 ds_y=ds_y,
                                                                 dtau=dtau,
                                                                 tau=0.0,
                                                                 s_x=x_lag,
                                                                 s_y=y_lag,
                                                                 horizontal_dist=False)

                    print(acf)
                    acfs[:,yi,xi]=acf
                    errs[:,yi,xi]=err
                    
            ho=h5py.File("mpi/%s/2d_mer_zon_res-%06d.h5"%(name,pi),"w")
            ho["acfs"]=acfs
            ho["errs"]=errs
            ho["s_x"]=s_x
            ho["s_y"]=s_y
            ho["ds_x"]=ds_x
            ho["ds_y"]=ds_y
            ho["dtau"]=dtau
            ho.close()
            
md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()


def summer():
    avg_mer_zon_acfs(md, 
                     dcos_thresh=0.8,
                     h0=90.0,
                     dh=2,
                     ds_x=50.0,
                     ds_y=50.0,
                     dtau=600.0,
                     s_x=n.arange(0.0, 400.0, 50.0),
                     s_y=n.arange(0.0, 400.0, 50.0),
                     ds_z=1,
                     years=[2018,2019,2020],
                     months=[5],
                     name="summer_2d_mer_zon_long")

summer()
