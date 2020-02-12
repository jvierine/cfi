#!/usr/bin/env python

import glob
import numpy as n
import matplotlib.pyplot as plt
import h5py
import datetime
import cfi_config as conf

fl=glob.glob("%s/*.h5"%(conf.data_directory))
fl.sort()
n_files=len(fl)
print(n_files)
links=["Alta_Alta","Andenes_Andenes","Andenes_Straumen","Tromso_Tromso","All"]
plot_all=False

for link in links:
    print(link)
    counts=[]
    countso=[]
    dates=[]
    
    for f in fl:
        h=h5py.File(f,"r")
        if link=="All":
            idx=n.arange(len(h["link"].value))
        else:
            idx=n.where(h["link"].value==link)[0]
        if len(idx)> 0:
            counts.append(len(n.unique(h["t"].value[idx])))
            print(n.unique(h["link"].value))
            countso.append(len(h["t"].value[idx]))
            dt_object = datetime.datetime.fromtimestamp(n.min(h["t"].value[idx]))
            dates.append(dt_object)
        h.close()
    if plot_all:
        plt.plot(dates,countso,label=link)
    else:
        plt.plot(dates,counts,label=link)

plt.legend()
plt.xlabel("Date (UTC)")
plt.ylabel("Counts per day")
if plot_all:
    plt.title("MMARIA Norway\ndetections per day")
else:
    plt.title("MMARIA Norway\nunique detections per day")

plt.show()
    

