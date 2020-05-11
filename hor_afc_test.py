import mmaria_read as mr
import cfi_dac as cfi
import cfi_config as c
import time
import datetime
import numpy as n
import h5py
import matplotlib.pyplot as plt
        

md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()
y=2019
y1=2019
m0=9
m1=12
day=1


    
d=md.read_data_date(d0=datetime.date(y,m0,day),d1=datetime.date(y1,m1,day))

meas=cfi.get_meas(meas_file=d,
             mean_rem=False,
             plot_dops=False,
             dcos_thresh=0.8,
             mean_wind_file="res/mean_wind_4h.h5",
             data='mmaria')


title='500km horisontal AFC with data from date {}.{}.{} to {}.{}.{}'.format(y,m0,day,y1,m1,day)
cfi.hor_acfs(meas,h0=90.0,dh=5,ds_z=1.0,ds_h=25.0,s_h=n.arange(0,500.0,50.0),dtau=600.0,title=title)

# In[] plotting made files

ho=h5py.File(title,"r")
shs=n.copy(ho["shs"])
acfs=n.copy(ho["acfs"])
    
dtau=n.copy(ho["dtau"])
ds_h=n.copy(ho["ds_h"])
ds_z=1.0

names=["$G_{uu}$","$G_{vv}$","$G_{ww}$",
       "$G_{uv}$","$G_{uw}$","$G_{vw}$"]

for i in range(6):
    plt.plot(shs,acfs[:,i],label=names[i])

plt.legend()
plt.xlabel("Horizontal lag (km)")    
plt.ylabel("Correlation (m$^2$/s$^2$)")
plt.title("Horizontal ACF $\Delta s_z=%1.1f$ km, $\Delta \\tau = %1.1f$ s $\Delta s_h=%1.1f$ km"%(ds_z,dtau,ds_h))
