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
                 name="summer_hacf",
                 remove_mean=False):

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
            pars.append([year,month,1])
    n_pars=len(pars)
    print("n_runs %d"%(n_pars))

    for pi in range(rank,n_pars,size):
        year=pars[pi][0]
        month=pars[pi][1]
        day=pars[pi][2]
        d0=datetime.date(year,month,day)
        t0=time.mktime(d0.timetuple())
        
        d=md.read_data(t0=t0,t1=t0+31*24*3600.0)
        n_meas=len(d["t"])
        print("n_meteors %d"%(n_meas))
        if n_meas > 100:
            if remove_mean:
                times,times_h,v,ve,rgs,lat0,lon0,dt,dh=mw.mean_wind(meas=d, 
                                                                    dt=mean_wind_time_avg,
                                                                    dh=1.0,
                                                                    max_alt=105,
                                                                    min_alt=78,
                                                                    dcos_thresh=dcos_thresh,
                                                                    ofname="res/tmp-%03d.h5"%(rank),
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
                                                      title=name)
            
            ho=h5py.File("mpi/%s/hacf_res-%06d.h5"%(name,pi),"w")
            ho["acfs"]=acfs
            ho["errs"]=errs
            ho["s_h"]=si_h
            ho["ds_h"]=ds_h
            ho["ds_z"]=ds_z
            ho["dtau"]=dtau
            ho.close()
            
# 2018 to 2020
# summer: 5,6,7
# winter: 11,12,1
# fall: 8,9,10
# spring: 2,3,4
#y=2019
#y1=2019
#m0=5
#m1=8
#day=1

# avg_hor_acfs(dcos_thresh=0.8,
#              mean_wind_time_avg=4*3600.0,
#              h0=90.0,
#              dh=5,
#              ds_h=25.0,
#              dtau=300.0,
#              s_h=n.arange(5.0,500.0,25.0),
#              ds_z=1,
#              years=[2018,2019,2020],
#              months=[5,6,7],
#              name="summer_hacf",
#              remove_mean=False,
#              n_days=31)

md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()

def winter():
    avg_hor_acfs(md,
                 dcos_thresh=0.8,
                 mean_wind_time_avg=4*3600.0,
                 h0=90.0,
                 dh=3,
                 ds_h=25.0,
                 dtau=600.0,
                 s_h=n.arange(0.0,500.0,25.0),
                 ds_z=1,
                 years=[2018,2019,2020],
                 months=[11,12,1],
                 name="winter_m_hacf",
                 remove_mean=False)
def winter_hp():
    avg_hor_acfs(md,
                 dcos_thresh=0.7,
                 mean_wind_time_avg=4*3600.0,
                 h0=90.0,
                 dh=2,
                 ds_h=25.0,
                 dtau=600.0,
                 s_h=n.arange(0.0,500.0,25.0),
                 ds_z=1,
                 years=[2018,2019,2020],
                 months=[11,12,1],
                 name="winter_hp_hacf",
                 remove_mean=True)
def spring():
    avg_hor_acfs(md,
                 dcos_thresh=0.8,
                 mean_wind_time_avg=4*3600.0,
                 h0=90.0,
                 dh=2,
                 ds_h=25.0,
                 dtau=600.0,
                 s_h=n.arange(0.0,500.0,25.0),
                 ds_z=1,
                 years=[2018,2019,2020],
                 months=[2,3,4],
                 name="spring_hacf",
                 remove_mean=False)
def fall():
    avg_hor_acfs(md,
                 dcos_thresh=0.8,
                 mean_wind_time_avg=4*3600.0,
                 h0=90.0,
                 dh=2,
                 ds_h=25.0,
                 dtau=600.0,
                 s_h=n.arange(0.0,500.0,25.0),
                 ds_z=1,
                 years=[2018,2019,2020],
                 months=[8,9,10],
                 name="fall_hacf",
                 remove_mean=False)

def summer():
    avg_hor_acfs(md,
                 dcos_thresh=0.8,
                 mean_wind_time_avg=4*3600.0,
                 h0=90.0,
                 dh=3,
                 ds_h=25.0,
                 dtau=600.0,
                 s_h=n.arange(0.0,500.0,25.0),
                 ds_z=1,
                 years=[2018,2019],
                 months=[5,6,7],
                 name="summer_m_hacf",
                 remove_mean=False)
def summer_hp():
    avg_hor_acfs(md,
                 dcos_thresh=0.7,
                 mean_wind_time_avg=4*3600.0,
                 h0=90.0,
                 dh=2,
                 ds_h=25.0,
                 dtau=600.0,
                 s_h=n.arange(0.0,500.0,25.0),
                 ds_z=1,
                 years=[2018,2019],
                 months=[5,6,7],
                 name="summer_hp_hacf",
                 remove_mean=True)

#exit(0)
#comm.Barrier()

#winter_hp()
#summer_hp()

#summer()
winter()
#spring()
#fall()
