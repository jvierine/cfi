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

data_prefix=""
if "data_prefix" in c["cfi"].keys():
    data_prefix=json.loads(c["cfi"]["data_prefix"])

# some parameters with default values
dcos_thresh=0.8
if "dcos_thresh" in c["cfi"].keys():
    dcos_thresh=float(json.loads(c["cfi"]["dcos_thresh"]))
    
debug_epsilon_fit=False
if "debug_epsilon_fit" in c["cfi"].keys():
    debug_epsilon_fit=bool(json.loads(c["cfi"]["debug_epsilon_fit"]))
    
high_pass_filter=True
if "high_pass_filter" in c["cfi"].keys():
    high_pass_filter=bool(json.loads(c["cfi"]["high_pass_filter"]))

horizontal_correlation_LTz=True
if "horizontal_correlation_LTz" in c["cfi"].keys():
    horizontal_correlation_LTz=bool(json.loads(c["cfi"]["horizontal_correlation_LTz"]))

debug_monthly_epsilon=False
if "debug_monthly_epsilon" in c["cfi"].keys():
    debug_monthly_epsilon=bool(json.loads(c["cfi"]["debug_monthly_epsilon"]))

debug_cfi_meteors=False
if "debug_cfi_meteors" in c["cfi"].keys():
    debug_cfi_meteors=bool(json.loads(c["cfi"]["debug_cfi_meteors"]))


horizontal_correlation_n_days=1
if "horizontal_correlation_n_days" in c["cfi"].keys():
    horizontal_correlation_n_days=float(json.loads(c["cfi"]["horizontal_correlation_n_days"]))

horizontal_correlation_post_avg=1
if "horizontal_correlation_post_avg" in c["cfi"].keys():
    horizontal_correlation_post_avg=int(json.loads(c["cfi"]["horizontal_correlation_post_avg"]))

epsilon_hlimit=[70,110]
if "epsilon_hlimit" in c["cfi"].keys():
    epsilon_hlimit=n.array(json.loads(c["cfi"]["epsilon_hlimit"]))
    
horizontal_correlation_heights=n.array(json.loads(c["cfi"]["horizontal_correlation_heights"]))
horizontal_correlation_dh=float(json.loads(c["cfi"]["horizontal_correlation_dh"]))

horizontal_correlation_dtau=600
if "horizontal_correlation_dtau" in c["cfi"].keys():
    horizontal_correlation_dtau=float(json.loads(c["cfi"]["horizontal_correlation_dtau"]))

horizontal_fit_synoptic=False
if "horizontal_fit_synoptic" in c["cfi"].keys():
    horizontal_fit_synoptic=bool(json.loads(c["cfi"]["horizontal_fit_synoptic"]))


horizontal_fit_min_lag=4
horizontal_fit_max_lag=18

if "horizontal_fit_max_lag" in c["cfi"].keys():
    horizontal_fit_max_lag=int(json.loads(c["cfi"]["horizontal_fit_max_lag"]))
if "horizontal_fit_min_lag" in c["cfi"].keys():
    horizontal_fit_min_lag=int(json.loads(c["cfi"]["horizontal_fit_min_lag"]))

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
