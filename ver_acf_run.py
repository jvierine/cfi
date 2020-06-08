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


def avg_ver_acfs(dcos_thresh=0.8,
                 mean_wind_time_avg=4*3600.0,
                 h0=80.0,
                 ds_h=100.0,
                 dtau=300.0,
                 s_z=n.arange(0.0,20.0,1.0),
                 years=[2018,2019,2020],
                 months=[5,6,7],
                 mean_rem=False,
                 name="summer_vacf",
                 n_days=31):

    os.system("mkdir -p mpi/%s"%(name))
    os.system("rm mpi/%s/*.h5"%(name))

    
    pars=[]
    for year in years:
        #for month in [5,6,7]:
        for month in months:
            pars.append((year,month,1))
    n_pars=len(pars)
    print("n_pars %d"%(n_pars))
    for pi in range(rank,n_pars,size):
        year=pars[pi][0]
        month=pars[pi][1]
        d0=datetime.date(year,month,1)
        t0=time.mktime(d0.timetuple())
    
        d=md.read_data(t0=t0,t1=t0+n_days*24*3600)
        n_meas=len(d["t"])
        print("n_meteors %d"%(n_meas))
        if n_meas > 100:
            mwname="res/temp-%f.h5"%(time.time())
            if mean_rem:
                times,times_h,v,ve,rgs,lat0,lon0,dt,dh=mw.mean_wind(meas=d, 
                                                                    dt=mean_wind_time_avg,
                                                                    dh=1.0,
                                                                    max_alt=105,
                                                                    min_alt=78,
                                                                    dcos_thresh=dcos_thresh,
                                                                    ofname=mwname,
                                                                    data='dict')
                
            meas=cfi.get_meas(meas_file=d,
                              mean_rem=mean_rem,
                              plot_dops=False,
                              dcos_thresh=dcos_thresh,
                              mean_wind_file=mwname,
                              data='mmaria')
            if mean_rem:
                os.system("rm %s"%(mwname))
            
            sz,acfs,errs,names,ds_z,tlag=cfi.ver_acfs(meas,
                                                      h0=h0,
                                                      dh=2.0,
                                                      s_z=s_z,
                                                      dtau=dtau,
                                                      tau=0.0,
                                                      s_h=0.0,
                                                      ds_h=ds_h)
                
            ho=h5py.File("mpi/%s/ver_res-%d-%f.h5"%(name,month,time.time()),"w")
            ho["acf"]=acfs
            ho["err"]=errs
            ho["s_z"]=s_z
            ho["ds_z"]=ds_z
            ho["ds_h"]=ds_h
            ho["dtau"]=dtau
            ho["h0"]=h0
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

#avg_ver_acfs(dcos_thresh=0.8,
#             mean_wind_time_avg=4*3600.0,
#             h0=80.0,
#             ds_h=100.0,
#             dtau=300.0,
#             s_z=n.arange(0.0,20.0,1.0),
#             years=[2018,2019,2020],
#             months=[5,6,7],
#             name="summer_vacf",
#             n_days=31)

def winter_large_scale():
    avg_ver_acfs(dcos_thresh=0.8,
                 mean_wind_time_avg=4*3600.0,
                 h0=88.0,
                 ds_h=500.0,
                 dtau=3600.0,
                 s_z=n.arange(0.0,20.0,1.0),
                 years=[2018,2019,2020],
                 months=[11,12,1],
                 name="winter_88_vacf_ls",
                 n_days=31)
    
def winter_small_scale():
    avg_ver_acfs(dcos_thresh=0.8,                 
                 h0=82.0,
                 ds_h=50.0,
                 dtau=600.0,
                 s_z=n.arange(0.0,20.0,1.0),
                 years=[2018,2019,2020],
                 months=[11,12,1],
                 name="winter_82_vacf",
                 n_days=31)
    
def summer_large_scale():
    avg_ver_acfs(dcos_thresh=0.8,
                 mean_wind_time_avg=4*3600.0,
                 h0=88.0,
                 ds_h=500.0,
                 dtau=3600.0,
                 s_z=n.arange(0.0,20.0,1.0),
                 years=[2018,2019,2020],
                 months=[11,12,1],
                 name="summer_88_vacf_ls",
                 n_days=31)
    
def summer_small_scale():
    avg_ver_acfs(dcos_thresh=0.8,
                 h0=82.0,
                 ds_h=50.0,
                 dtau=600.0,
                 s_z=n.arange(0.0,20.0,1.0),
                 years=[2018,2019,2020],
                 months=[5,6,7],
                 name="summer_82_vacf",
                 n_days=31)

summer_small_scale()
#summer_large_scale()
#winter_small_scale()
#winter_large_scale()

