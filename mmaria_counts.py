#!/usr/bin/env python

import glob
import numpy as n
import matplotlib.pyplot as plt
import h5py
import datetime

fl=glob.glob("*.h5")
fl.sort()
n_files=len(fl)
print(n_files)
links=[b"Alta_Alta",b"Andenes_Andenes",b"Andenes_Straumen",b"Tromso_Tromso"]
plot_all=False

for link in links:
    print(link)
    counts=[]
    countso=[]
    dates=[]
    
    for f in fl:
        h=h5py.File(f,"r")

        idx=n.where(h["link"][()]==link)[0]
        if len(idx)> 0:
            counts.append(len(n.unique(h["t"][()][idx])))
            countso.append(len(h["t"][()][idx]))
            dt_object = datetime.datetime.fromtimestamp(n.min(h["t"][()][idx]))
            dates.append(dt_object)

        h.close()
        
    if plot_all:
        plt.plot(dates,countso,label=str(link))
    else:
        plt.plot(dates,counts,label=str(link))

plt.legend()
plt.xlabel("Date (UTC)")
plt.ylabel("Counts per day")
if plot_all:
    plt.title("MMARIA Norway\ndetections per day")
else:
    plt.title("MMARIA Norway\nunique detections per day")

plt.show()
    

