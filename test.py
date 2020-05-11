
import mmaria_read as mr
import cfi_dac as cfi
import cfi_config as c
import time
import datetime
import numpy as n

md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()

d=md.read_data_yyyymmdd(d0=datetime.date(2019,1,5),d1=datetime.date(2019,2,1))

meas=cfi.get_meas(meas_file=d,
             mean_rem=False,
             plot_dops=False,
             dcos_thresh=0.8,
             mean_wind_file="res/mean_wind_4h.h5",
             data='mmaria')



dtau=1800.0
h_max=1*24.0
n_t=int(h_max*3600.0/dtau)
    
    
cfi.temporal_acfs(meas,h0=93.5,dh=2,ds_z=1.0,ds_h=25.0,dtau=dtau,tau=n.arange(float(n_t))*dtau)


