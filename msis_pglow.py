import pyglow
import scipy.constants as sc

from datetime import datetime

masses={"HE":4*sc.atomic_mass,"O":16*sc.atomic_mass,"N2":28*sc.atomic_mass,"O2":32*sc.atomic_mass,"AR":40*sc.atomic_mass,"H":sc.atomic_mass,"N":14*sc.atomic_mass,"O_anomalous":0}


def get_atmospheric_mass(year,month,hour,alt=80.0,lat=69.0,lon=19.0):
    dn=datetime(year,month,1,hour,0,1)

    pt=pyglow.Point(dn,lat,lon,alt)
    pt.run_msis()

    constituents = pt.nn.keys()
#    print(pt.nn)
    mass=0.0
    for c in constituents:
        # mass kg/m^2
        mass+=pt.nn[c]*masses[c]*1e6
    return(mass)

if __name__ == "__main__":
    print(get_atmospheric_mass(2010,6,hour=23,alt=90))
