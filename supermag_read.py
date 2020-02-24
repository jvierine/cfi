#!/usr/bin/env python
#
# read supermag style csv data
#
import numpy as n
import matplotlib.pyplot as plt
import datetime
import time

import supermag_stations as ss

class supermag_data:
    def __init__(self,fname="test_bl.csv",plot=False):
        sms = ss.supermag_stations()
        
        f=open(fname,"r")
        stats=[]
        dtimes=[]
        neus=[]
        lats=[]
        lons=[]
        for li,l in enumerate(f.readlines()):
            if li > 0:
                line=l.split(",")
                date=line[0]
                d=datetime.datetime.strptime("%sZ"%(date), "%Y-%m-%d %H:%M:%SZ")
                unixtime = time.mktime(d.timetuple())
                dtimes.append(unixtime)
                #        print(d)
                #       print(unixtime)
                stat=line[1]

                p=sms.prop(stat)
                lats.append(p["glat"])
                lons.append(p["glon"])
                stats.append(stat)
                neu=n.array([float(line[6]),float(line[7]),float(line[8])])
                neus.append(neu)
                #        print(ned)
                #       print(li)
                #      print(l)
                #stat_names=n.unique(stats)
        dtimes=n.array(dtimes)
        neus=n.array(neus)
        stats=n.array(stats)
        self.lons=n.array(lons)
        self.lats=n.array(lats)
        if plot:
            idx=n.where(stats == 'ABK')[0]
            plt.plot(dtimes[idx],neus[idx,0])
            plt.plot(dtimes[idx],neus[idx,1])
            plt.plot(dtimes[idx],neus[idx,2])
            plt.show()
            print(idx)
        self.times=dtimes
        self.neu=neus
        self.stats=stats

        
    def get_bounds(self):
        return([n.min(self.times),n.max(self.times)])

    
    def get_meas(self,t0,t1,station=None):
        if station == None:
            idx=n.where((self.times >= t0)&(self.times <= t1))[0]
        else:
            idx=n.where((self.times >= t0)&(self.times <= t1)&(self.stats==station))[0]
        return({"times":self.times[idx],
                "neu":self.neu[idx,:],
                "stat":self.stats[idx],
                "glat":self.lats[idx],
                "glon":self.lons[idx]})
    
#    def get_meas_station(self,t0,t1,station="TRO"):
 #       m=self.get_meas(t0,t1)
  #      idx=n.where(m["stat"]==station)[0]
   #     return({"times":times[idx],
    #            "neu":times[idx],
     #           "stat":times[idx],
      #          "glat":times[idx],
       #         "glon":times[idx],

if __name__ == "__main__":
    d=supermag_data(fname="tro_2018_2019.csv",plot=False)
    b=d.get_bounds()
    # get all data from "TRO"
    m=d.get_meas(b[0],b[1],station="TRO")
    plt.plot((m["times"]-b[0])/3600.0,m["neu"][:,0])
    plt.show()
    
