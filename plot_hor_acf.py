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

def estimate_epsilons():
    n_heights=len(c.horizontal_correlation_heights)
    
    
    E=n.zeros([12,n_heights])

    hidx=n.argsort(c.horizontal_correlation_heights)
    heights=c.horizontal_correlation_heights[hidx]

    rhidx=[]
    for h in c.horizontal_correlation_heights:
        rhidx.append(n.argmin(n.abs(heights - h)))
    rhidx=n.array(rhidx,dtype=n.int)
    
  #  print(hidx)
 #   print(heights)
#    exit(0)
    for hi,h0 in enumerate(c.horizontal_correlation_heights):
        
        name="hcor_%1.2f"%(h0)#94.00"
        print("calculating %s"%(name))

        fl=glob.glob("%s/%s/*.h5"%(c.data_directory,name))
        fl.sort()
        if len(fl) == 0:
            # if no results yet, then skip
            continue

#        print(fl)
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

        if c.debug_monthly_epsilon:
            print(all_acfs.shape)
            plt.subplot(121)
            days=n.arange(all_acfs.shape[0])
            plt.pcolormesh(s_h,t0s,all_acfs[:,:,0],vmin=0,vmax=500,cmap="jet")
#            plt.title("%1.0f km"%(h))
            plt.colorbar()
            plt.xlabel("Horizontal separation (km)")
            plt.ylabel("Date")
            plt.title("$R'_{LL}$ %1.0f km"%(h0))
            plt.subplot(122)
            plt.pcolormesh(s_h,t0s,all_acfs[:,:,1],vmin=0,vmax=500,cmap="jet")
            plt.colorbar()
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

            if average > 0:
                for j in range(-average,average):
                    if ( (j+i) >= 0) and ( (j+i) < all_acfs.shape[0]):
                        n_avg+=1.0
                        avg_Ruu+=all_acfs[j+i,2:21,0]
                        avg_Rvv+=all_acfs[j+i,2:21,1]
            else:
                n_avg=1.0
                avg_Ruu+=all_acfs[i,2:21,0]
                avg_Rvv+=all_acfs[i,2:21,1]

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

        if c.debug_monthly_epsilon:
            plt.plot(t0s,eps_uu*1e3,".")
            plt.ylim([-10,100])
            plt.xlabel("Date")
            plt.ylabel("Dissipation rate (mW/kg)")
            plt.tight_layout()
            #plt.plot(t0s,eps_vv*1e3,".")
            #plt.plot(t0s,0.5*(eps_uu+eps_vv)*1e3,".")
            plt.show()
            
        # reverse lookup index to array of heights
        E[:,rhidx[hi]]=eps_month/eps_count
            
        if c.debug_monthly_epsilon:
            plt.plot(n.arange(1,13),eps_month/eps_count,".-")
            plt.ylim([0,50])
            plt.xlabel("Month")
            plt.ylabel("Dissipation rate (mW/kg)")
            plt.tight_layout()
            plt.show()

    
    plt.pcolormesh(n.arange(12),heights,n.transpose(E),vmin=10,vmax=30)
    plt.colorbar()
    plt.show()

    
if __name__ == "__main__":
    estimate_epsilons()
