#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import mmaria_read as mr
import cfi_config as c

def get_height_histogram(n_t=None,d_delta=7):
    """
    Look at the height histograms.
    """
    md=mr.mmaria_data(c.data_directory)

    b=md.get_bounds()

    if n_t == None:
        n_t=int(n.floor((b[1]-b[0])/(d_delta*24*3600)))
    print(n_t)


    t0=24*3600*n.floor(b[0]/(24*3600))

    n_ranges=30
    rgs=n.linspace(70,110,num=n_ranges+1)

    S=n.zeros([n_t,n_ranges])

    for di in range(n_t):
        d=md.read_data(t0+di*d_delta*24*3600,t0+(di+1)*d_delta*24*3600)
        hist,bins=n.histogram(d["heights"]/1e3,bins=rgs)
        ranges=0.5*(bins[0:len(hist)]+bins[1:(len(hist)+1)])
        S[di,:]=hist

    tvec=d_delta*n.arange(n_t)        
    plt.pcolormesh(tvec,ranges,n.transpose(S))
    plt.plot(tvec,ranges[n.argmax(S,axis=1)])
    plt.show()


if __name__ == "__main__":
    get_height_histogram(d_delta=7)
