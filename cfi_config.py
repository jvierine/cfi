#!/usr/bin/env python

from configparser import ConfigParser
import sys
import json

if len(sys.argv) != 2:
    print("provide config file argument")
    
c = ConfigParser()
c.read(sys.argv[1])

lat=float(json.loads(c["cfi"]["lat"]))
lon=float(json.loads(c["cfi"]["lon"]))

data_directory=json.loads(c["cfi"]["data_directory"])
plot_directory=json.loads(c["cfi"]["plot_directory"])

#data_directory="/data0/SIMONe_multilink/JRO"
#plot_directory="/data0/SIMONe_multilink/JRO/"

#lat=c["cfi"]["lat"]
#data_directory="D:\Data\mmaria_data"

# latitude and longitude of the middle point of the network
# only used for map centering in plot, so it is not critical
# that this is correct.
#lat=0.0
#lon=19.0


