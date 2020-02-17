#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import mmaria_read as mr
import cfi_config as c

def get_decay_histogram(n_t=None,d_delta=1,height=[70,110],logalpha=[0,2.0],N=50):
    """
    Look at the height histograms.
    """
    md=mr.mmaria_data(c.data_directory)

    b=md.get_bounds()

    if n_t == None:
        n_t=int(n.floor((b[1]-b[0])/(d_delta*24*3600)))
        
    hbins=n.linspace(height[0],height[1],num=N)
    abins=n.linspace(logalpha[0],logalpha[1],num=N)
    S=n.zeros([N,N])
    
    t0=24*3600*n.floor(b[0]/(24*3600))

    for di in range(n_t):
        d=md.read_data(t0+di*d_delta*24*3600,t0+(di+1)*d_delta*24*3600)
        bragg_scale=n.sqrt(d["braggs"][:,0]**2.0+d["braggs"][:,1]**2.0+d["braggs"][:,2]**2.0)
        print(bragg_scale)
        H,x,y=n.histogram2d(n.log10(d["alpha_norm"]*bragg_scale**2.0),d["heights"]/1e3,bins=[abins,hbins])        
        plt.pcolormesh(abins,hbins,n.transpose(H))
        plt.show()

if __name__ == "__main__":
    get_decay_histogram(d_delta=31)
