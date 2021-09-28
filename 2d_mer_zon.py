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

def avg_mer_zon(md, # data
                dcos_thresh=0.8,
                h0=87.0,         # smallest height
                h1=93.0,         # biggest height
                ds_x=25.0,       # longest separation in east-west dist km
                ds_y=25.0,       # longest separation in north-south dist km
                ds_z=1,          # longest separation in vertical km 
                dtau=1800.0,     # longest separation in time (s)
                s_x=0.0,         # x lag mean
                s_y=0.0,         # y lag mean 
                xi=0,            # index 
                yi=0,            # index
                x_lags=n.arange(-500,500.0,50.0),
                y_lags=n.arange(-500,500.0,50.0),
                years=[2018,2019,2020],                 
                months=[5,6,7],
                n_months=3,
                name="summer_2d_mer_zon"):

    print("zon %1.2f mer %1.2f"%(s_x,s_y))

    dh=(93.0-87.0)*2.0

    # to be done. divide into arbitrary time intervals
    pars=[]
    for year in years:
        for month in months:
            pars.append([year,month,1])
    
    n_pars=len(pars)
    print("n_runs %d"%(n_pars))

    for pi in range(n_pars):
        year=pars[pi][0]
        month=pars[pi][1]
        day=pars[pi][2]
        d0=datetime.date(year,month,day)
        t0=time.mktime(d0.timetuple())        
        d=md.read_data(t0=t0,t1=t0+n_months*31*24*3600.0,h0=h0,h1=h1)
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
                                                         h0=0.5*(h0+h1),
                                                         dh=dh,
                                                         ds_z=1.0,
                                                         ds_x=ds_x,
                                                         ds_y=ds_y,
                                                         dtau=dtau,
                                                         tau=0.0,
                                                         s_x=s_x,
                                                         s_y=s_y,
                                                         horizontal_dist=False)
            print(ixs)
            print(iys)
            print(acf)
            ofname="mpi/%s/2d_mer_zon_res-%d-%d-%d-%d.h5"%(name,xi,yi,pi,rank)
            print("writing %s"%(ofname))
            ho=h5py.File(ofname,"w")
            ho["acf"]=acf
            ho["err"]=err
            ho["s_x"]=x_lags
            ho["s_y"]=y_lags
            ho["x"]=s_x
            ho["y"]=s_y
            ho["xi"]=xi
            ho["yi"]=yi
            ho["ds_x"]=ds_x
            ho["ds_y"]=ds_y
            ho["dtau"]=dtau
            ho.close()
            
md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()

def mpi_run(name="all_2d"):

    h0=87
    h1=93
    s_x=n.arange(-350, 350.0, 50.0)
    s_y=n.arange(-350, 350.0, 50.0)
    ds_x = 50.0
    ds_y = 50.0
    dtau = 1800.0

    if rank == 0:
        os.system("rm mpi/%s/*.h5"%(name))        
        os.system("mkdir -p mpi/%s"%(name))
    comm.Barrier()

    pars=[]
    for yi,ylag in enumerate(s_y):
        for xi,xlag in enumerate(s_x):
            pars.append((yi,xi,ylag,xlag))
    n_par=len(pars)
    print(n_par)

    for pi in range(rank,n_par,size):
        p=pars[pi]
        yi=p[0]
        xi=p[1]
        ylag=p[2]
        xlag=p[3]
        print(p)
        avg_mer_zon(md, 
                    dcos_thresh=0.8,
                    h0=h0,
                    h1=h1,                    
                    ds_x=ds_x,
                    ds_y=ds_y,
                    dtau=dtau,
                    s_y=ylag,
                    s_x=xlag,
                    y_lags=s_y,
                    x_lags=s_x,
                    yi=yi,
                    xi=xi,
                    ds_z=1,
                    years=[2018,2019],
                    months=[1,2,3,4,5,6,7,8,9,10,11,12],
                    n_months=1,
                    name=name)

    print("done %d"%(rank))

#mpi_run(name="summer_mer_zon_100km_900s",month0=5)
mpi_run(name="all")
#mpi_run(name="winter_mer_zon_50km_900s",month0=11)
#mpi_run(name="winter_mer_zon_50km_900s",month0=11)
