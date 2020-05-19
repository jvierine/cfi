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
                      h0=90.0,
                      dh=2,
                      ds_h=50.0,
                      dtau=300.0,
                      taus=n.arange(0.0,7*24*3600,300.0),
                      tau=0.0,
                      ti=0,
                      ds_z=1,
                      years=[2018,2019,2020],
                      months=[5],
                      n_days=31*3,
                      name="summer_tacf"):
    if rank == 0:
        os.system("rm mpi/%s/*.h5"%(name))
        os.system("mkdir -p mpi/%s"%(name))
    comm.Barrier()

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
    
        d=md.read_data(t0=t0,t1=t0+n_days*24*3600.0,h0=h0-dh,h1=h0+dh)
        n_meas=len(d["t"])
        print("n_meteors %d"%(n_meas))
        if n_meas > 100:
            meas=cfi.get_meas(meas_file=d,
                              mean_rem=False,
                              plot_dops=False,
                              dcos_thresh=dcos_thresh,
                              data='dict')

            
            acf, err, itaus, ixs, iys, izs, is_h=cfi.cfi(meas,
                                                         h0=h0,
                                                         dh=dh,
                                                         ds_z=1.0,
                                                         ds_h=ds_h,
                                                         s_h=0.0,
                                                         dtau=dtau,
                                                         tau=tau,
                                                         horizontal_dist=True)
            print(acf)
            ho=h5py.File("mpi/%s/tacf_res-%03d-%06d.h5"%(name,ti,pi),"w")
            ho["acf"]=acf
            ho["err"]=err
            ho["ti"]=ti
            ho["tau"]=tau
            ho["taus"]=taus
            ho["dtau"]=dtau
            ho["ds_h"]=ds_h
            ho["h0"]=h0
            ho.close()

def summer(name="summer_tacf_14_50km_300",h0=90):
    
    if rank == 0:
        os.system("rm mpi/%s/*.h5"%(name))        
        os.system("mkdir -p mpi/%s"%(name))
    comm.Barrier()

    tau=n.arange(0,14*24*3600,300.0)
    pars=[]
    for ti,tlag in enumerate(tau):
        pars.append((ti,tlag))
        
    n_par=len(pars)
    print(n_par)

    for pi in range(rank,n_par,size):
        p=pars[pi]
        ti=p[0]
        tlag=p[1]
    
        avg_temporal_acfs(dcos_thresh=0.8,
                          h0=h0,
                          dh=2,
                          ds_h=50.0,
                          dtau=300.0,
                          taus=tau,
                          tau=tlag,
                          ti=ti,
                          years=[2018,2019,2020],
                          months=[5],
                          n_days=31*3,
                          name=name)
    

#winter_small_scale()
summer()
