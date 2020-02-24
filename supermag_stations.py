#!/usr/bin/env python
#
# Gjerloev, J. W. (2012), The SuperMAG data processing technique, J. Geophys. Res., 117, A09213, doi:10.1029/2012JA017683.
# Gjerloev, J. W. (2009), A Global Ground-Based Magnetometer Initiative, EOS, 90, 230-231, doi:10.1029/2009EO270002.
#
import numpy as n
import matplotlib.pyplot as plt

class supermag_stations():
    def __init__(self,fname="20200221-09-06-supermag-stations.csv"):
        
        self.props={}
        f=file(fname,"r")
        for l in f.readlines():
            s=l.split(",")
            name=s[0]
            if name != "IAGA":
                prop = {"name":s[0],
                        "glon":float(s[1]),
                        "glat":float(s[2]),
                        "mlon":float(s[3]),
                        "mlat":float(s[4]),
                        "name":s[5],
                        "net":s[7],
                        "opnum":s[6]}
                self.props[name]=prop
        f.close()
        
    def prop(self,name):
        return(self.props[name])

if __name__ == "__main__":
    s=supermag_stations()
    print(s.prop("TRO"))
