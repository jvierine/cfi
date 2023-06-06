#!/usr/bin/env python

from configparser import ConfigParser
import sys
import json
import numpy as n

if len(sys.argv) != 2:
    print("provide config file argument")
    
c = ConfigParser()
c.read(sys.argv[1])

lat=float(json.loads(c["cfi"]["lat"]))
lon=float(json.loads(c["cfi"]["lon"]))

data_directory=json.loads(c["cfi"]["data_directory"])
plot_directory=json.loads(c["cfi"]["plot_directory"])

# some parameters with default values
dcos_thresh=0.8
if "dcos_thresh" in c["cfi"].keys():
    dcos_thresh=float(json.loads(c["cfi"]["dcos_thresh"]))
    
debug_epsilon_fit=False
if "debug_epsilon_fit" in c["cfi"].keys():
    debug_epsilon_fit=bool(json.loads(c["cfi"]["debug_epsilon_fit"]))

horizontal_correlation_n_days=1
if "horizontal_correlation_n_days" in c["cfi"].keys():
    horizontal_correlation_n_days=float(json.loads(c["cfi"]["horizontal_correlation_n_days"]))

horizontal_correlation_post_avg=1
if "horizontal_correlation_post_avg" in c["cfi"].keys():
    horizontal_correlation_post_avg=int(json.loads(c["cfi"]["horizontal_correlation_post_avg"]))
    
horizontal_correlation_heights=n.array(json.loads(c["cfi"]["horizontal_correlation_heights"]))
horizontal_correlation_dh=float(json.loads(c["cfi"]["horizontal_correlation_dh"]))


#data_directory="/data0/SIMONe_multilink/JRO"
#plot_directory="/data0/SIMONe_multilink/JRO/"

#lat=c["cfi"]["lat"]
#data_directory="D:\Data\mmaria_data"

# latitude and longitude of the middle point of the network
# only used for map centering in plot, so it is not critical
# that this is correct.
#lat=0.0
#lon=19.0


if __name__ == "__main__":
    print(c)
