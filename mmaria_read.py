#!/usr/bin/env python

import numpy as n
import glob
import h5py
import matplotlib.pyplot as plt
import cfi_config as c
import time
import datetime


#alpha_norm               Dataset {56452}
#braggs                   Dataset {56452, 3}
#dcos                     Dataset {56452, 2}
#dh                       Dataset {SCALAR}
#dop_errs                 Dataset {56452}
#dops                     Dataset {56452}
#dt                       Dataset {SCALAR}
#heights                  Dataset {56452}
#lats                     Dataset {56452}
#link                     Dataset {56452}
#lons                     Dataset {56452}
#rgs                      Dataset {30}
#t                        Dataset {56452}
#times                    Dataset {48}
#v                        Dataset {2, 48, 30}
#v_resid                  Dataset {56452}
#ve                       Dataset {2, 48, 30}

keys=["alpha_norm","braggs","dcos","dh","dop_errs","dops",
      "dt","heights","lats","link","lons","rgs","t","rgs",
      "t","times","v","v_resid","ve"]


class mmaria_data:
    def __init__(self,dname,debug=False):
        self.fl = glob.glob("%s/*.h5"%(dname))
        self.fl.sort()
        self.mint=[]
        self.maxt=[]
        self.debug=debug
        
        for f in self.fl:
            if self.debug:
                print("reading %s"%(f))
            h=h5py.File(f,"r")
            t=h["t"].value
            self.mint.append(n.min(t))
            self.maxt.append(n.max(t))
            h.close()
        self.mint=n.array(self.mint)
        self.maxt=n.array(self.maxt)        

    def get_bounds(self):
        return([n.min(self.mint),n.max(self.maxt)])

    def read_data_date(self,d0,d1,read_all_detections=True):#is datetime.date(d0) valid to run a function, so can just put in 2019,1,1
        #d0,d1 is a date like d0=datetime.date(2019,1,5)
        
        t0=time.mktime(d0.timetuple())
        t1=time.mktime(d1.timetuple())
        return(self.read_data(t0,t1,read_all_detections))
        
    def read_data(self,
                  t0,
                  t1,
                  h0=0,
                  h1=150,
                  read_all_detections=True):
        """
        Read all meteor radar network data between these times (unix)
        """
        file_idx=[]
        for fi in range(len(self.fl)):
            if self.maxt[fi] > t0 and self.mint[fi] < t1:
                file_idx.append(fi)

        alpha_norm=n.zeros([0],dtype=n.float32)
        braggs=n.zeros([0,3],dtype=n.float32)
        dcos=n.zeros([0,2],dtype=n.float32)
        dh=1.5
        dop_errs=n.zeros([0],dtype=n.float32)
        dops=n.zeros([0],dtype=n.float32)
        dt=0
        heights=n.zeros([0],dtype=n.float32)
        lats=n.zeros([0],dtype=n.float32)
        lons=n.zeros([0],dtype=n.float32)        
        link=n.zeros([0],dtype="<U60")
        rgs=n.zeros([30],dtype=n.float32)
        t=n.zeros([0],dtype=n.float32)
        times=n.zeros([0],dtype=n.float32)
        v=n.zeros([2,0,30],dtype=n.float32)
        v_resid=n.zeros([0],dtype=n.float32)
        ve=n.zeros([2,0,30],dtype=n.float32)
        print("reading")
        for i in file_idx:#range(first_idx,last_idx+1):
#            print(i)
            h=h5py.File(self.fl[i],"r")
            if read_all_detections:
                didx=n.where( ((h["t"].value) > t0) &
                              ((h["t"].value) < t1) &
                              (h["heights"].value/1e3 > h0) &
                              (h["heights"].value/1e3 < h1) )[0]
                
                t = n.concatenate((t,h["t"].value[didx]))
                alpha_norm = n.concatenate((alpha_norm,h["alpha_norm"].value[didx]))
                braggs = n.concatenate((braggs,h["braggs"].value[didx,:]))
        
                dcos = n.concatenate((dcos,h["dcos"].value[didx,:]))
                dh=h["dh"].value
                dop_errs = n.concatenate((dop_errs,h["dop_errs"].value[didx]))
                dops = n.concatenate((dops,h["dops"].value[didx]))
                dt=h["dt"].value
                heights=n.concatenate((heights,h["heights"].value[didx]))
                lats=n.concatenate((lats,h["lats"].value[didx]))
                lons=n.concatenate((lons,h["lons"].value[didx]))
                link=n.concatenate((link,h["link"].value[didx]))
                v_resid=n.concatenate((v_resid,h["v_resid"].value[didx]))
            
            # mean horizontal wind model
            rgs=h["rgs"].value
            times=n.concatenate((times,h["times"].value))
            v=n.concatenate((v,h["v"].value),axis=1)
            ve=n.concatenate((ve,h["ve"].value),axis=1)
            
            if self.debug:
                print("file idx %d"%(i))
        tu,idx=n.unique(t,return_index=True)
        return({"t":t[idx],
                "alpha_norm":alpha_norm[idx],
                "braggs":braggs[idx,:],
                "dcos":dcos[idx,:],
                "dh":dh,
                "dop_errs":dop_errs[idx],
                "dops":dops[idx],
                "dt":dt,
                "heights":heights[idx],
                "lats":lats[idx],
                "lons":lons[idx],
                "link":link[idx],
                "rgs":rgs,
                "v_resid":v_resid,
                "times":times,
                "v":v,
                "ve":ve})


if __name__ == "__main__":
    """
    Example usage
    """
    # directory with all mmaria network data
    md=mmaria_data(c.data_directory)
    # what is the data bounds (first and last time stamp)
    print(md.get_bounds())
    
    # read all meteor radar data between these two timestamps
    d=md.read_data(1514774804,1514974804)

    plt.pcolormesh(d["times"],d["rgs"]/1e3,n.transpose(d["v"][0,:,:]),vmin=-100,vmax=100)
    plt.xlabel("Time (unix)")
    plt.ylabel("Altitude (km)")
    plt.colorbar()
    plt.show()

