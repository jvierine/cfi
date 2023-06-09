import time
import datetime
import numpy as n
import h5py
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

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
    ED=n.zeros([365,n_heights])
    ED[:,:]=n.nan
    EE=n.zeros([365,n_heights])
    EE[:,:]=n.nan

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
        t00=[]        
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
                t00.append(h["t0"][()])
                n_avg+=1.0
            else:
                print("too many nans in acf. aborting")
            h.close()


        all_acfs=n.array(all_acfs)

        if c.debug_monthly_epsilon:
            print(all_acfs.shape)
            fig=plt.figure(layout="constrained")
            gs=GridSpec(nrows=3,
                        ncols=2,
                        width_ratios=[1,1],
                        height_ratios=[0.05,1,1],
                        figure=fig)
            gs.update(left=0.05,right=0.95,bottom=0.08,top=0.93,wspace=0.02,hspace=0.03)
            ax1=fig.add_subplot(gs[1,0])
            ax2=fig.add_subplot(gs[0:2,1])
            ax3=fig.add_subplot(gs[2,0:2])
            cbax=fig.add_subplot(gs[0,0])
#            plt.subplot(221)
            days=n.arange(all_acfs.shape[0])
            im=ax1.pcolormesh(s_h[1:],t0s,all_acfs[:,1:,0],vmin=0,vmax=500,cmap="jet")
            cb=plt.colorbar(im,cax=cbax,orientation="horizontal",ticklocation="top")
            cb.set_label("$R_{LL}$ (m$^2$/s$^2$)")
            #ax1.colorbar()
#            plt.colorbar(p,cax=ax1)
            ax1.set_xlabel("$s$ (km)")
            ax1.set_ylabel("Date")
#            ax1.set_title("$R_{LL}$")
#            plt.title("%1.0f km"%(h))
#            plt.colorbar()
 #           plt.xlabel("Horizontal separation (km)")
  #          plt.ylabel("Date")
   #         plt.title("$R'_{LL}$ %1.0f km"%(h0))
#            plt.subplot(222)
#            p=ax2.pcolormesh(s_h,t0s,all_acfs[:,:,1],vmin=0,vmax=500,cmap="jet")
 #           plt.colorbar(p,cax=ax2)
 #           ax2.set_xlabel("Horizontal separation (km)")
  #          ax2.set_ylabel("Date")
   #         ax2.set_title("$R'_{TT}$")
    #        plt.colorbar()
     #       plt.xlabel("Horizontal separation (km)")
      #      plt.ylabel("Date")
       #     plt.title("$R'_{TT}$")
#            plt.tight_layout()
 #           plt.show()
 
#            fig.suptitle("%1.0f km"%(h0))

        # average this many days before and after
        # essentially a one month rolling avg
        average = c.horizontal_correlation_post_avg

        eps_uu=[]
        eps_vv=[]
        eps_uu_std=[]
        eps_vv_std=[]        

        eps_month=n.zeros(12)
        eps_count=n.zeros(12)
        t0_doy=n.zeros(all_acfs.shape[0])

        min_lag=4
        max_lag=18
        n_lags=max_lag-min_lag
        
        resids_uu=n.zeros([all_acfs.shape[0],n_lags])
        resids_vv=n.zeros([all_acfs.shape[0],n_lags])        
        for i in range(all_acfs.shape[0]):

            avg_Ruu=n.zeros(n_lags)
            avg_Rvv=n.zeros(n_lags)

            n_avg=0.0

            if average > 0:
                for j in range(-average,average):
                    if ( (j+i) >= 0) and ( (j+i) < all_acfs.shape[0]):
                        n_avg+=1.0
                        avg_Ruu+=all_acfs[j+i,min_lag:max_lag,0]
                        avg_Rvv+=all_acfs[j+i,min_lag:max_lag,1]
            else:
                n_avg=1.0
                avg_Ruu+=all_acfs[i,min_lag:max_lag,0]
                avg_Rvv+=all_acfs[i,min_lag:max_lag,1]

            avg_Ruu=avg_Ruu/n_avg
            avg_Rvv=avg_Rvv/n_avg    

            fr=fit_epsilon.fit_epsilon0(s_h[min_lag:max_lag],avg_Ruu,debug=c.debug_epsilon_fit)

            # plot equinox and summer
            summer_idx=87
            equi_idx=154
            if (i==summer_idx) or (i==equi_idx) :
                if c.debug_monthly_epsilon:
                    color="C0"
                    if i==summer_idx:
                        color="C1"
                    ax2.plot(s_h[min_lag:max_lag],fr["model"],label="$\epsilon=%1.0f\pm%1.1f$ mW/kg"%(fr["xhat"][0]*1e3,2*fr["eps_std"]*1e3),color=color)
                    ax2.errorbar(s_h[min_lag:max_lag],avg_Ruu,yerr=2*fr["err_std"],fmt=".",color=color)
                    ax2.grid()
                    ax2.set_ylim([0,350])
                    #ax2.set_title("
#                    ax2.set_title("%1.0f km"%(h0))
               #     ax2.set_title("$\epsilon=%1.2f \pm %1.2f$ mW/kg"%(fr["xhat"][0]*1e3,fr["eps_std"]*1e3))
                    #        plt.plot(s,R,".")
                    ax2.set_xlabel("$s$ (km)")
                    ax2.set_ylabel("$R_{LL}$ (m$^2$/s$^2$)")

                    
                #p=ax2.pcolormesh(s_h,t0s,all_acfs[:,:,1],vmin=0,vmax=500,cmap="jet")

            
            xhat_uu=fr["xhat"]
            res_uu=fr["resid"]
            eps_uu_std.append(fr["eps_std"])
            resids_uu[i,:]=res_uu
            fr=fit_epsilon.fit_epsilon0(s_h[min_lag:max_lag],avg_Rvv,debug=c.debug_epsilon_fit)
            xhat_vv=fr["xhat"]
            res_vv=fr["resid"]
            resids_vv[i,:]=res_vv
            
            eps_uu.append(xhat_uu[0])
            eps_vv.append(xhat_vv[0])

            eps_month[int(t0s[i].month-1)]+=xhat_uu[0]*1e3
            eps_count[int(t0s[i].month-1)]+=1.0

            
        eps_uu=n.array(eps_uu)
        eps_uu_std=n.array(eps_uu_std)        
        eps_vv=n.array(eps_vv)

        ED[n.array(n.floor(365*(t00-n.min(t00))/(24*3600*365)),dtype=n.int),rhidx[hi]]=eps_uu*1e3
        EE[n.array(n.floor(365*(t00-n.min(t00))/(24*3600*365)),dtype=n.int),rhidx[hi]]=eps_uu_std*1e3
        
        if c.debug_monthly_epsilon:
            ax2.legend()
#            plt.subplot(223)
            t00=n.array(t00)
            sidx=n.arange(0,len(t00),7)
            ax3.errorbar( 12*(t00[sidx]-n.min(t00))/(24*3600*365) ,eps_uu[sidx]*1e3, yerr=2*eps_uu_std[sidx]*1e3 ,fmt=".",color="black")
            ax3.set_ylabel("$\\varepsilon$ (mW/kg)")
            ax3.set_xlabel("Month")

#            plt.ylim([-10,100])
 #           plt.xlabel("Date")
  #          plt.ylabel("Dissipation rate (mW/kg)")
#            plt.tight_layout()
            #plt.plot(t0s,eps_vv*1e3,".")
            #plt.plot(t0s,0.5*(eps_uu+eps_vv)*1e3,".")
 #           plt.show()
            
        # reverse lookup index to array of heights
        E[:,rhidx[hi]]=eps_month/eps_count
            
        if c.debug_monthly_epsilon:
#            plt.subplot(224)
            ax3.plot(n.arange(12)+0.5,eps_month/eps_count,".-",color="red")
            ax3.grid()
            #         plt.ylim([0,50])
            #        plt.xlabel("Month")
            #       plt.ylabel("Dissipation rate (mW/kg)")
            plt.tight_layout()
            plt.show()

        if c.debug_monthly_epsilon:
            plt.subplot(221)
            plt.pcolormesh(resids_uu)
            plt.subplot(222)
            plt.pcolormesh(resids_vv)
            plt.subplot(223)
            plt.plot(n.sqrt(n.nanmean(n.abs(resids_uu)**2.0,axis=0)))
            plt.subplot(224)
            plt.plot(n.sqrt(n.nanmean(n.abs(resids_vv)**2.0,axis=0)))
            plt.show()

#    plt.subplot(121)
 #   plt.pcolormesh(12*n.arange(365)/365,heights,n.transpose(EE),vmin=0,vmax=50)
  #  plt.xlabel("Month")
 #   plt.ylabel("Height (km)")
  #  cb=plt.colorbar()
   # cb.set_label("$\\varepsilon$ error std (mW/kg)")
    
  #  plt.subplot(122)
    # don't show epsilon with error larger than 10 mW/kg
    ED[EE > 10]=n.nan
    months=12*n.arange(365)/365
    #heights
    mm,hh=n.meshgrid(months,heights)
    plt.pcolormesh(mm,hh,n.transpose(ED),vmin=0,vmax=50,cmap="jet")
    plt.ylim([79,97])
    plt.xlabel("Month")
    plt.ylabel("Height (km)")
    cb=plt.colorbar()
    cb.set_label("$\\varepsilon$ (mW/kg)")
    plt.tight_layout()
    plt.show()

    
if __name__ == "__main__":
    estimate_epsilons()
