import time
import datetime
import numpy as n
import h5py
import matplotlib.pyplot as plt

# our internal modules
import mmaria_read as mr
import cfi_dac as cfi
import cfi_config as c
import mean_wind_est as mw
import glob

def read_acf(name="summer_88_vacf"):
    
    fl=glob.glob("mpi/%s/*.h5"%(name))
    fl.sort()
    
    print(fl)
    h=h5py.File(fl[0],"r")
    acf=n.copy(h["acf"].value)
    err=n.copy(h["err"].value)
    s_z=n.copy(h["s_z"].value)
    dtau=n.copy(h["dtau"].value)
    h0=n.copy(h["h0"].value)
    #ds_h=n.copy(h["ds_h"].value)
    h.close()
    acf[:,:]=0
    err[:,:]=0
    all_acfs=[]
    all_errs=[]
    all_sz=[]
    n_avg=0
    for f in fl:
        h=h5py.File(f,"r")
        #   plt.plot(h["acf"].value[:,0])
        #    plt.show()
        if n.sum(n.isnan(h["acf"].value)) == 0:
            all_acfs.append(n.copy(h["acf"].value))
            all_errs.append(n.copy(h["err"].value))
            all_sz.append(n.copy(h["s_z"].value))
            n_avg+=1
        h.close()
    all_acfs=n.array(all_acfs)
    all_errs=n.array(all_errs)
    all_sz=n.array(all_sz)
    err_vars=n.zeros([len(s_z),6])
    acfs=n.zeros([len(s_z),6])
    s_z[:]=0.0
    for i in range(len(s_z)):
        s_z[i]=n.mean(all_sz[:,i])
        for ci in range(6):
            err_vars[i,ci]=n.var(all_acfs[:,i,ci])
            ws=0.0
            for mi in range(int(n_avg)):
                if n.sum(n.isnan(all_acfs[mi,i,ci]))== 0:
                    w=1.0/all_errs[mi,i,ci]
                    ws+=w
                    acfs[i,ci]+=w*all_acfs[mi,i,ci]
                else:
                    print("nans")
            acfs[i,ci]=(1.0/ws)*acfs[i,ci]
            err_vars[i,ci]=1/ws

    return(s_z,acfs,err_vars)

s_z,acfs_h,errs= read_acf(name="summer_82_vacf")
#s_z,acfs_l,errs_l= read_acf(name="summer_vacf_ls")

acfs=acfs_h#-acfs_l
plt.subplot(121)    
plt.plot(s_z,acfs_h[:,0],label=cfi.cf_names[0])
plt.plot(s_z,acfs_h[:,1],label=cfi.cf_names[1])
#plt.plot(s_z,acfs_h[:,1])
#plt.plot(s_z,acfs_l[:,1])
plt.xlabel("Vertical lag (km)")
plt.ylabel("Correlation function (m$^2$/s$^2$)")
plt.title("Autocorrelation function")
plt.subplot(122)    
plt.loglog(s_z,2*1.2*acfs[1,0]-2*acfs[:,0])
plt.loglog(s_z,2*1.2*acfs[1,1]-2*acfs[:,1])
plt.loglog(s_z,100.0*s_z**(2.0/3.0))
plt.xlabel("Vertical lag (km)")
plt.ylabel("Structure function (m$^2$/s$^2$)")
plt.title("Structure function")
plt.show()
