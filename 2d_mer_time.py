#
# 2d correlation function in meridional or zonal direction.
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
                 ti=0,
                 yi=0,
                 s_y=0.0,
                 s_x=0.0,
                 tau=0.0,
                 ds_z=1,
                 years=[2018,2019,2020],
                 y_lags=[],
                 t_lags=[],
                 month=5,
                 n_months=3,
                 name="summer_2d_mer_time"):


    
    pars=[]
    for year in years:
        pars.append([year,month,1])
    
    n_pars=len(pars)
    print("n_runs %d"%(n_pars))

    for pi in range(n_pars):
        year=pars[pi][0]
        month=pars[pi][1]
        day=pars[pi][2]
        d0=datetime.date(year,month,day)
        t0=time.mktime(d0.timetuple())
        
        d=md.read_data(t0=t0,t1=t0+n_months*31*24*3600.0, h0=h0-dh, h1=h0+dh)
        n_meas=len(d["t"])
        print("n_meteors %d"%(n_meas))
        
        if n_meas > 100:
            meas=cfi.get_meas(meas_file=d,
                              mean_rem=False,
                              plot_dops=False,
                              dcos_thresh=dcos_thresh,
                              mean_wind_file="res/tmp-%03d.h5"%(rank),
                              data='mmaria')
            

            acf, err, itaus, ixs, iys, izs, is_h=cfi.cfi(meas,
                                                         h0=h0,
                                                         dh=2,
                                                         ds_z=1.0,
                                                         ds_x=ds_x,
                                                         ds_y=ds_y,
                                                         dtau=dtau,
                                                         tau=tau,
                                                         s_x=s_x,
                                                         s_y=s_y,
                                                         horizontal_dist=False)
            print("y %1.2f t %1.2f"%(s_y,tau))
            print(acf)
            ho=h5py.File("mpi/%s/2d_mer_time_res-%03d-%03d-%06d.h5"%(name,yi,ti,pi),"w")
            ho["acfs"]=acf
            ho["errs"]=err
            ho["s_y"]=s_y
            ho["s_x"]=s_x
            ho["tau"]=tau
            ho["y_lags"]=y_lags
            ho["t_lags"]=t_lags
            ho["ti"]=ti
            ho["yi"]=yi
            ho["ds_x"]=ds_x
            ho["ds_y"]=ds_y
            ho["dtau"]=dtau
            ho.close()
            
md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()


def mpi_run(mer=True,name="summer_2d_zon_time_3",
            month0=5,
            dx=50.0):
    
    if rank == 0:
        os.system("rm mpi/%s/*.h5"%(name))        
        os.system("mkdir -p mpi/%s"%(name))
    comm.Barrier()

    s_y=n.arange(0,400.0, 50.0)
    tau=n.arange(0,12*3600,900.0)
    pars=[]
    for yi,ylag in enumerate(s_y):
        for ti,tlag in enumerate(tau):
            pars.append((yi,ti,ylag,tlag))
    n_par=len(pars)
    print(n_par)

    for pi in range(rank,n_par,size):
        p=pars[pi]
        yi=p[0]
        ti=p[1]
        ylag=p[2]
        tlag=p[3]

        if mer:
            sy_lag=ylag
            sx_lag=0.0
        else:
            sy_lag=0.0
            sx_lag=ylag
        
        avg_mer_acfs(md, 
                     dcos_thresh=0.8,
                     h0=90.0,
                     dh=2,
                     ds_x=dx,
                     ds_y=dx,
                     dtau=900.0,
                     s_y=sy_lag,
                     s_x=sx_lag,
                     tau=tlag,
                     y_lags=s_y,
                     t_lags=tau,
                     yi=yi,
                     ti=ti,
                     ds_z=1,
                     years=[2018,2019,2020],
                     month=month0,
                     n_months=3,
                     name=name)

#mpi_run(mer=True,name="winter_2d_mer_time_12_900s_50km",month0=11)
#mpi_run(mer=False,name="winter_2d_zon_time_12_900s_50km",month0=11)
mpi_run(mer=False,name="winter_2d_zon_time_12_900s_100km",month0=11,dx=100.0)
#mpi_run(mer=True,name="summer_2d_mer_time_12_900s_50km",month0=5)
#mpi_run(mer=False,name="summer_2d_zon_time_12_900s_50km",month0=5)
