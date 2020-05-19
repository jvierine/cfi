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

name="summer_hacf"
#name="winter_hacf"
#name="summer_hp_hacf"

fl=glob.glob("mpi/%s/*.h5"%(name))
fl.sort()

print(fl)
h=h5py.File(fl[0],"r")
acf=n.copy(h["acfs"].value)
err=n.copy(h["errs"].value)
s_h=n.copy(h["s_h"].value)
dtau=n.copy(h["dtau"].value)
ds_z=n.copy(h["ds_z"].value)
ds_h=n.copy(h["ds_h"].value)
h.close()
acf[:,:]=0
err[:,:]=0
all_acfs=[]
all_errs=[]
all_sh=[]
n_avg=0.0
for f in fl:
    h=h5py.File(f,"r")
 #   plt.plot(h["acfs"].value[:,0])
#    plt.show()
    if n.sum(n.isnan(h["acfs"].value)) == 0:
        all_acfs.append(n.copy(h["acfs"].value))
        all_errs.append(n.copy(h["errs"].value))
        all_sh.append(n.copy(h["s_h"].value))
        n_avg+=1.0
    h.close()
    
    
all_acfs=n.array(all_acfs)
all_errs=n.array(all_errs)
all_sh=n.array(all_sh)
err_vars=n.zeros([len(s_h),6])
acfs=n.zeros([len(s_h),6])
s_h[:]=0.0
for i in range(len(s_h)):
    s_h[i]=n.mean(all_sh[:,i])
    for ci in range(6):
        err_vars[i,ci]=n.var(all_acfs[:,i,ci])
        ws=0.0
        for mi in range(int(n_avg)):
            if n.sum(n.isnan(all_acfs[mi,i,ci]))== 0:

                w=1.0/all_errs[mi,i,ci]
#                print(w)
                ws+=w
                acfs[i,ci]+=w*all_acfs[mi,i,ci]
            else:
                print("nans")
#        print(ws)
        acfs[i,ci]=(1.0/ws)*acfs[i,ci]

                    
names=["$G_{uu}$","$G_{vv}$","$G_{ww}$",
       "$G_{uv}$","$G_{uw}$","$G_{vw}$"]
colors=["C0","C1","C2","C3","C4","C5"]
cfi.plot_hor_acfs(shs=s_h,
                  names=names,
                  acfs=acfs,
                  ds_z=1.0,
                  dtau=dtau,
                  ds_h=ds_h,
                  err_vars=err_vars,
                  colors=colors,
                  zlag=1.1,
                  n_avg=n_avg)

    # plt.savefig("C:/Users/OleK/Master_thesis/figs/fig_%s.png"%(name))
    # ho=h5py.File("res/%s.h5"%(name),"w")
    # ho["acf"]=acfs
    # ho["s_h"]=s_h
    # ho["err_var"]=err_vars/n.sqrt(n_avg)
    # ho["h0"]=h0
    # ho["dtau"]=dtau
    # ho["ds_h"]=ds_h
    # ho.close()
    
    
    
    
