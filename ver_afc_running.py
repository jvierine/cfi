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
                 name="summer_vacf",
                 n_days=31):
    n_lags=len(s_z)
    all_acfs=[]
    
    n_avg=0.0
    #for year in [2020,2018,2019]:
    for year in years:
        #for month in [5,6,7]:
        for month in months:
            d0=datetime.date(year,month,1)
            t0=time.mktime(d0.timetuple())
    
            for day in range(n_days):
                d=md.read_data(t0=t0+day*24*3600.0,t1=t0+(day+1)*24*3600.0)
                n_meas=len(d["t"])
                print("n_meteors %d"%(n_meas))
                if n_meas > 100:
    
                    times,times_h,v,ve,rgs,lat0,lon0,dt,dh=mw.mean_wind(meas=d, 
                                                                        dt=mean_wind_time_avg,
                                                                        dh=1.0,
                                                                        max_alt=105,
                                                                        min_alt=78,
                                                                        dcos_thresh=dcos_thresh,
                                                                        ofname="res/tmp.h5",
                                                                        data='dict')
                
                    meas=cfi.get_meas(meas_file=d,
                                      mean_rem=True,
                                      plot_dops=False,
                                      dcos_thresh=dcos_thresh,
                                      mean_wind_file="res/tmp.h5",
                                      data='mmaria')
            
                    #title='Vertical AFC with data from date {}.{}.{} to {}.{}.{}'.format(y,m0,day,y1,m1,day)
                    sz,acfs,errs,names,ds_z,tlag=cfi.ver_acfs(meas,h0=h0,dh=2.0,s_z=s_z,
                                                              dtau=dtau, tau=0.0, s_h=0.0, 
                                                              ds_h=ds_h)
                    all_acfs.append(acfs)#+=acfs
                    n_avg+=1.0
    all_acfs=n.array(all_acfs)
    err_vars=n.zeros([len(s_z),6])
    acfs=n.zeros([len(s_z),6])
    for i in range(len(s_z)):
        for ci in range(6):
            err_vars[i,ci]=n.var(all_acfs[:,i,ci])
            acfs[i,ci]=n.mean(all_acfs[:,i,ci])
                    
    #print(all_acfs.shape)
    colors=["C0","C1","C2","C3","C4","C5"]
    
    cfi.plot_ver_acf(szs=s_z,acfs=acfs,names=names,ds_z=ds_z,
                     dtau=dtau,err_vars=err_vars,colors=colors,
                     n_avg=n_avg)
    

    plt.show()
    plt.savefig("C:/Users/OleK/Master_thesis/figs/fig_%s.png"%(name))
    ho=h5py.File("res/%s.h5"%(name),"w")
    ho["acf"]=acfs
    ho["s_z"]=s_z
    ho["err_var"]=err_vars/n.sqrt(n_avg)
    ho["h0"]=h0
    ho["dtau"]=dtau
    ho["ds_h"]=ds_h
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

avg_ver_acfs(dcos_thresh=0.8,
             mean_wind_time_avg=4*3600.0,
             h0=80.0,
             ds_h=100.0,
             dtau=300.0,
             s_z=n.arange(0.0,20.0,1.0),
             years=[2018,2019,2020],
             months=[11,12,1],
             name="winter_vacf",
             n_days=31)
