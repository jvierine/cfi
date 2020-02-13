#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import mmaria_read as mr
import cfi_config as c

def get_horizontal_histogram(n_t=None,d_delta=1,lat_min=65,lat_max=74,lon_min=5,lon_max=32,N=50):
    """
    Look at the height histograms.
    """
    md=mr.mmaria_data(c.data_directory)

    b=md.get_bounds()

    if n_t == None:
        n_t=int(n.floor((b[1]-b[0])/(d_delta*24*3600)))
        
    lons=n.linspace(lon_min,lon_max,num=N)
    lats=n.linspace(lat_min,lat_max,num=N)
    S=n.zeros([N,N])
    
    t0=24*3600*n.floor(b[0]/(24*3600))

    for di in range(n_t):
        d=md.read_data(t0+di*d_delta*24*3600,t0+(di+1)*d_delta*24*3600)
        H,x,y=n.histogram2d(d["lons"],d["lats"],bins=[lons,lats])
        plt.pcolormesh(lons,lats,n.transpose(H))
        plt.show()


if __name__ == "__main__":
    get_horizontal_histogram(d_delta=60)
