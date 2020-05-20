from mpl_toolkits.basemap import Basemap

import numpy as n
import matplotlib.pyplot as plt
import cfi_config as c
import mmaria_read as mr
import stuffr


#works for mmaria_data
md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()

t0=stuffr.date2unix(2019,1,1,0,0,0)
t1=stuffr.date2unix(2019,1,2,0,0,0)

# about two weeks of data in June 2019
d=md.read_data(t0,t1)

# direction cosine filter
dcos_thresh=0.8
lats=d["lats"]
lons=d["lons"]
dcoss=d["dcos"]
dcos2=n.sqrt(dcoss[:,0]**2.0+dcoss[:,1]**2.0)
ok_idx=n.where(dcos2 < dcos_thresh)[0]

lat0=n.median(lats)
lon0=n.median(lons)

m=Basemap(projection="stere",
          lat_0=lat0,
          lon_0=lon0,
          llcrnrlat=66,
          urcrnrlat=72,
          llcrnrlon=11,
          urcrnrlon=32,
          resolution="h")
m.drawcoastlines()

pars=n.arange(67,72,2)
m.drawparallels(pars,labels=pars)

mers=n.arange(15,30,5)
m.drawmeridians(mers,labels=mers)

x,y=m(lons[ok_idx],lats[ok_idx])
plt.plot(x,y,".",alpha=0.1)
plt.show()
