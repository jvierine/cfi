import numpy as n
import cfi_config as c

Re=6371.0
latdeg2km=n.sin(n.pi/180.0)*Re
londeg2km=n.pi*Re*n.cos(n.pi*c.lat/180.0)/180.0#65.122785 
