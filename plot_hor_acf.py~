import time
import datetime
import numpy as n
import h5py
import matplotlib.pyplot as plt

import fit_epsilon

# our internal modules
import mmaria_read as mr
import cfi_dac as cfi
import cfi_config as c
import mean_wind_est as mw
import glob
import stuffr

# 600 s, 90 km 4 hour wind removed
#name="daily"

# 300 s, 90 km 4 hour wind removed
#name="daily_300s_3day_90km"

#name="daily_600s_5day_90km_LTz"
#name="daily_600s_5day_90km_LTz"
#name="daily_600s_5day_82km_LTz"

#name="peru_daily_600s_7day_90km_LTz"


name="hcor_90.00"

fl=glob.glob("%s/%s/*.h5"%(c.data_directory,name))
fl.sort()

print(fl)
h=h5py.File(fl[0],"r")
acf=n.copy(h["acfs"][()])
err=n.copy(h["errs"][()])
s_h=n.copy(h["s_h"][()])
dtau=n.copy(h["dtau"][()])
ds_z=n.copy(h["ds_z"][()])
ds_h=n.copy(h["ds_h"][()])
h.close()
acf[:,:]=0
err[:,:]=0
all_acfs=[]
all_errs=[]
all_sh=[]
n_avg=0.0


#months = n.arange(12)+1

#all_acfs=n.zeros([12,])
t0s=[]
for f in fl:
    h=h5py.File(f,"r")
    print(h["acfs"][()])
 #   plt.plot(h["acfs"][()][:,0])
#    plt.show()
#
 #   if n.sum(n.isnan(h["acfs"][()])) == 0:
    if True:
        all_acfs.append(n.copy(h["acfs"][()]))
        all_errs.append(n.copy(h["errs"][()]))
        all_sh.append(n.copy(h["s_h"][()]))
        t0s.append(stuffr.unix2date(h["t0"][()]))
        n_avg+=1.0
    else:
        print("too many nans in acf. aborting")
    h.close()
    
    
all_acfs=n.array(all_acfs)
print(all_acfs.shape)
plt.subplot(121)
days=n.arange(all_acfs.shape[0])
plt.pcolormesh(s_h,t0s,all_acfs[:,:,0],vmin=0,vmax=400)
plt.xlabel("Horizontal separation (km)")
plt.ylabel("Date")
plt.title("$R'_{LL}$")
plt.subplot(122)
plt.pcolormesh(s_h,t0s,all_acfs[:,:,1],vmin=0,vmax=400)
plt.xlabel("Horizontal separation (km)")
plt.ylabel("Date")
plt.title("$R'_{TT}$")
plt.tight_layout()
plt.show()


# average this many days before and after
# essentially a one month rolling avg
average = c.horizontal_correlation_post_avg

eps_uu=[]
eps_vv=[]

eps_month=n.zeros(12)
eps_count=n.zeros(12)
t0_doy=n.zeros(all_acfs.shape[0])
for i in range(all_acfs.shape[0]):

    avg_Ruu=n.zeros(19)
    avg_Rvv=n.zeros(19)

    n_avg=0.0

    for j in range(-average,average):
        if ( (j+i) >= 0) and ( (j+i) < all_acfs.shape[0]):
            n_avg+=1.0
            avg_Ruu+=all_acfs[j+i,2:21,0]
            avg_Rvv+=all_acfs[j+i,2:21,1]

    avg_Ruu=avg_Ruu/n_avg
    avg_Rvv=avg_Rvv/n_avg    
    
    xhat_uu=fit_epsilon.fit_epsilon0(s_h[2:21],avg_Ruu,debug=c.debug_epsilon_fit)
    xhat_vv=fit_epsilon.fit_epsilon0(s_h[2:21],avg_Rvv,debug=c.debug_epsilon_fit)
    eps_uu.append(xhat_uu[0])
    eps_vv.append(xhat_vv[0])

    eps_month[int(t0s[i].month-1)]+=xhat_uu[0]*1e3
    eps_count[int(t0s[i].month-1)]+=1.0
    
eps_uu=n.array(eps_uu)
eps_vv=n.array(eps_vv)
plt.plot(t0s,eps_uu*1e3,".")
plt.ylim([0,100])
plt.xlabel("Date")
plt.ylabel("Dissipation rate (mW/kg)")
plt.tight_layout()
#plt.plot(t0s,eps_vv*1e3,".")
#plt.plot(t0s,0.5*(eps_uu+eps_vv)*1e3,".")
plt.show()

plt.plot(n.arange(1,13),eps_month/eps_count,".-")
plt.ylim([0,50])
plt.xlabel("Month")
plt.ylabel("Dissipation rate (mW/kg)")
plt.tight_layout()

plt.show()

#exit(0)
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
                ws+=w
                acfs[i,ci]+=w*all_acfs[mi,i,ci]
            else:
                print("nans")
        acfs[i,ci]=(1.0/ws)*acfs[i,ci]

# todo: fit the zero lag using an exp function



        
#print(acfs)
#plt.loglog(s_h,2.0*acfs[1,0] - 2*acfs[:,0])
#plt.show()

# fix zero lag!
#acfs[0,:]=acfs[1,:]
                    
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
    
    
    
    
