import glob
import numpy as n
import matplotlib.pyplot as plt
import h5py
import datetime

import cfi_config as c
import mmaria_read as mr



def meas_time_of_day(meas, hour_of_day=n.arange(48)*0.5,
                     dtau=1800.0):
    
    time=meas["t"]

    n_times=int((n.max(time)-n.min(time))/(dtau))
    t0=n.min(time)
    #y-axis=number of meas during the 30 min interval
    #x-axis=30 min time bins
    
    hod=n.mod(time/3600.0,24)
    # histogram and interpolate the number of measurements as a function of day
    hods=n.array(hod)
    # 30 minute bins, 0..24 hours utc histogram 
    thist, tbins=n.histogram(n.mod(hods/3600.0,24),bins=48)
    plt.figure(0)
    plt.bar(hour_of_day, thist)
    title='Time of day when meteor measurements occur'
    plt.title(title)
    plt.xlabel('hour of day')
    plt.ylabel('number of measurments')
    plt.savefig("C:/Users/OleK/Master_thesis/figs/fig_%s.png"%(title))
    plt.xticks(n.arange(0,25, step=2))
    plt.show()

# In[] making plot for # meas on altitude

def meas_altitude(meas, links):
    for link in links:
        print(link)
    
        if link == "All":
            idx=n.arange(len(d["link"]),dtype=n.int)
        else:
            idx=n.where(d["link"]==link)[0]
        if len(idx)> 0:
            print(n.unique(d["link"]))
            
            heights=n.around(d["heights"][idx],-2)
            alts, counts_alts = n.unique(heights, return_counts=True)
            
            plt.figure(1)                   #plots # of measerments on height
            plt.plot(counts_alts, alts, label=link)
            
    title='Number of measurements on altitude'    
    plt.xlabel("# of measurements")
    plt.ylabel("Altitude")
    plt.title(title)
    plt.legend()
    plt.savefig("C:/Users/OleK/Master_thesis/figs/fig_%s.png"%(title))



md=mr.mmaria_data(c.data_directory)#for many files in a directory
b=md.get_bounds()

d=md.read_data_date(d0=datetime.date(2019,1,1),d1=datetime.date(2019,2,1))

links=["Alta_Alta","Andenes_Andenes","Andenes_Straumen","Tromso_Tromso","All"]

#meas_time_of_day(meas=d)#makes bar plot for # of meas on time of day
meas_altitude(meas=d,links=links)