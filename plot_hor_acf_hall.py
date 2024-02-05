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

import os

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
    EDvv=n.zeros([365,n_heights])
    EDvv[:,:]=n.nan
    EE=n.zeros([365,n_heights])
    EE[:,:]=n.nan
    S=n.zeros([365,n_heights])
    S[:,:]=n.nan
    R0S=n.zeros([365,n_heights])
    R0S[:,:]=n.nan
    R0Svv=n.zeros([365,n_heights])
    R0Svv[:,:]=n.nan

    hidx=n.argsort(c.horizontal_correlation_heights)
    heights=c.horizontal_correlation_heights[hidx]

    rhidx=[]
    for h in c.horizontal_correlation_heights:
        rhidx.append(n.argmin(n.abs(heights - h)))
    rhidx=n.array(rhidx,dtype=int)
    
  #  print(hidx)
 #   print(heights)
#    exit(0)
    for hi,h0 in enumerate(c.horizontal_correlation_heights):

        if c.debug_monthly_epsilon:                    
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

        
        name="hcor_%1.2f"%(h0)#94.00"
        print("calculating %s"%(name))

        fl=glob.glob("%s/%s/%s/*.h5"%(c.data_directory,c.data_prefix,name))
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
        all_aacfs=n.copy(all_acfs)


        # average this many days before and after
        # essentially a one month rolling avg
        average = c.horizontal_correlation_post_avg

        eps_uu=[]
        ke_uu=[]
        ke_vv=[]                
        eps_vv=[]
        eps_uu_std=[]
        eps_vv_std=[]
        synop_uu=[]
        synop_vv=[]                        

        eps_month=n.zeros(12)
        eps_count=n.zeros(12)
        t0_doy=n.zeros(all_acfs.shape[0])

        min_lag=c.horizontal_fit_min_lag#4
        max_lag=c.horizontal_fit_max_lag#18
        n_lags=max_lag-min_lag
        
        resids_uu=n.zeros([all_acfs.shape[0],n_lags])
        resids_vv=n.zeros([all_acfs.shape[0],n_lags])
        
        for i in range(all_acfs.shape[0]):

            avg_Ruu=n.zeros(all_acfs.shape[1])
            avg_Rvv=n.zeros(all_acfs.shape[1])

            n_avg=0.0

            if average > 0:
                for j in range(-average,average):
                    if ( (j+i) >= 0) and ( (j+i) < all_acfs.shape[0]):
                        n_avg+=1.0
                        avg_Ruu+=all_acfs[j+i,:,0]
                        avg_Rvv+=all_acfs[j+i,:,1]
            else:
                n_avg=1.0
                avg_Ruu+=all_acfs[i,:,0]
                avg_Rvv+=all_acfs[i,:,1]

            avg_Ruu=avg_Ruu/n_avg
            avg_Rvv=avg_Rvv/n_avg

            all_aacfs[i,:,0]=avg_Ruu
            all_aacfs[i,:,1]=avg_Rvv

            
            fr=fit_epsilon.fit_epsilon0(s_h[min_lag:max_lag],avg_Ruu[min_lag:max_lag],debug=c.debug_epsilon_fit,fit_synoptic=c.horizontal_fit_synoptic)

            if c.debug_monthly_epsilon:            
                # plot equinox and summer
                summer_idx=91#87
                equi_idx=182#154
                if (i==summer_idx) or (i==equi_idx) :
                    if c.debug_monthly_epsilon:
                        color="C0"
                        if i==summer_idx:
                            color="C1"
#                        ax2.plot(s_h[min_lag:max_lag],fr["model"],label="$\epsilon=%1.0f$ mW/kg DoY=%d"%(fr["xhat"][0]*1e3,i),color=color)
                        ax2.plot(s_h[min_lag:max_lag],fr["model"],label="DoY=%d"%(i),color=color)                        
                        ax2.errorbar(s_h[min_lag:max_lag],avg_Ruu[min_lag:max_lag],yerr=2*fr["err_std"],fmt=".",color=color)
                        ax2.grid()
                        ax2.set_ylim([10,180])
                        #ax2.set_title("
                        ax2.set_title("%1.0f km"%(h0))
                   #     ax2.set_title("$\epsilon=%1.2f \pm %1.2f$ mW/kg"%(fr["xhat"][0]*1e3,fr["eps_std"]*1e3))
                        #        plt.plot(s,R,".")
                        ax2.set_xlabel("$s$ (km)")
                        ax2.set_ylabel("$R_{L}$ (m$^2$/s$^2$)")

                    
                #p=ax2.pcolormesh(s_h,t0s,all_acfs[:,:,1],vmin=0,vmax=500,cmap="jet")

            
            xhat_uu=fr["xhat"]
            res_uu=fr["resid"]
            eps_uu_std.append(fr["eps_std"])
            resids_uu[i,:]=res_uu
            fr=fit_epsilon.fit_epsilon0(s_h[min_lag:max_lag],avg_Rvv[min_lag:max_lag],debug=c.debug_epsilon_fit,fit_synoptic=c.horizontal_fit_synoptic)
            xhat_vv=fr["xhat"]
            res_vv=fr["resid"]
            resids_vv[i,:]=res_vv

            # if the zero-lag estimate is more than the std
            # otherwise consider the result to be too noisy
            #            if fr["xhat"][1] > 2*fr["err_std"][0]:
            # avg of the L and T dissipation rates.
            eps_uu.append(xhat_uu[0])
            eps_vv.append(xhat_vv[0])
            ke_uu.append(xhat_uu[1])
            ke_vv.append(xhat_vv[1])

            if c.horizontal_fit_synoptic:
                synop_uu.append(xhat_uu[2])
                synop_vv.append(xhat_vv[2])                
            #else:
             #   eps_uu.append(n.nan)
              #  eps_vv.append(n.nan)
               # ke_uu.append(n.nan)
                #ke_vv.append(n.nan)                

            eps_month[int(t0s[i].month-1)]+=xhat_uu[0]*1e3
            eps_count[int(t0s[i].month-1)]+=1.0



        if c.debug_monthly_epsilon:
            print(all_acfs.shape)
#            plt.subplot(221)



            days=n.arange(all_acfs.shape[0])
            use_colormesh=True
            if use_colormesh:
                im=ax1.pcolormesh(s_h[1:],t0s,all_aacfs[:,1:,0],vmin=0,vmax=400)
                cb=plt.colorbar(im,cax=cbax,orientation="horizontal",ticklocation="top")
                cb.set_label("$R_{L}$ (m$^2$/s$^2$)")
                #ax1.colorbar()
                #            plt.colorbar(p,cax=ax1)
                ax1.set_xlabel("$s$ (km)")
                ax1.set_ylabel("Date")
            else:
 #               for di in range(all_acfs.shape[0]):
#                    ax1.plot(s_h[:],all_acfs[di,:,0],".")
                ax1.plot(s_h,n.mean(all_acfs[:,:,0],axis=0))
                ax1.plot(s_h,n.mean(all_acfs[:,:,1],axis=0))                
                    
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

            
        eps_uu=n.array(eps_uu)
        eps_uu_std=n.array(eps_uu_std)        
        eps_vv=n.array(eps_vv)        
        ke_uu=n.array(ke_uu)
        ke_vv=n.array(ke_vv)

        synop_uu=n.array(synop_uu)
        synop_vv=n.array(synop_vv)        

        
        ED[n.array(n.floor(365*(t00-n.min(t00))/(24*3600*365)),dtype=int),rhidx[hi]]=eps_uu*1e3
        EDvv[n.array(n.floor(365*(t00-n.min(t00))/(24*3600*365)),dtype=int),rhidx[hi]]=eps_vv*1e3        
        
        R0S[n.array(n.floor(365*(t00-n.min(t00))/(24*3600*365)),dtype=int),rhidx[hi]]=ke_uu
        R0Svv[n.array(n.floor(365*(t00-n.min(t00))/(24*3600*365)),dtype=int),rhidx[hi]]=ke_vv        
        EE[n.array(n.floor(365*(t00-n.min(t00))/(24*3600*365)),dtype=int),rhidx[hi]]=eps_uu_std*1e3

        if c.horizontal_fit_synoptic:        
            S[n.array(n.floor(365*(t00-n.min(t00))/(24*3600*365)),dtype=int),rhidx[hi]]=10**synop_uu
        
        if c.debug_monthly_epsilon:
            ax2.legend()
#            plt.subplot(223)
            t00=n.array(t00)
#            sidx=n.arange(0,len(t00),int(c.horizontal_correlation_n_days+c.horizontal_correlation_post_avg+1))
            sidx=n.arange(0,len(t00),int(c.horizontal_correlation_n_days))
            ax3.errorbar( 365*(12*(t00[sidx]-n.min(t00))/(24*3600*365))/12.0 ,eps_uu[sidx]*1e3, yerr=2*eps_uu_std[sidx]*1e3 ,fmt=".",color="black",label="3-day avg")
            if os.path.exists("data/hallmf.h5"):
                hh=h5py.File("data/hallmf.h5","r")
                #plt.plot(hh["doy"][()],hh["dr"][()]/4.0,".")
                ax3.plot(hh["doy"][()],hh["dr"][()]/4.0,label="$0.25 \cdot \\varepsilon_{\mathrm{MF}}$ (Hall et.al., 99)",color="gray",zorder=2)
                hh.close()
            ax3.set_ylabel("$\\varepsilon$ (mW/kg)")
            ax3.set_xlabel("Day of year")


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
            ax3.plot(365*(n.arange(12)+0.5)/12.0,eps_month/eps_count,".-",color="red",zorder=3,label="Monthly avg")
            ax3.legend(ncol=3)            
            ax3.grid()
            ax3.set_ylim([0,35])
            #         plt.ylim([0,50])
            #        plt.xlabel("Month")
            #       plt.ylabel("Dissipation rate (mW/kg)")
            plt.tight_layout()
            plt.savefig("figs/hprof-%d.png"%(h0),dpi=400)
            plt.show()

            ho=h5py.File("data/eps-%1.2f.h5"%(h0),"w")
            ho["eps"]=eps_uu
            ho["eps_std"]=eps_uu_std
            ho["ke_uu"]=ke_uu
            ho["t"]=t00
            ho.close()
            

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

    # average values of epsilon and sigma after fitting for them
    # this causes less bias than averaging horizontal correlation functions
    if c.epsilon_post_avg > 1:
        w=n.repeat(1.0/c.epsilon_post_avg,c.epsilon_post_avg)
        for rgi in range(ED.shape[1]):
            ED[:,rgi]=n.convolve(ED[:,rgi],w,mode="same")
            R0S[:,rgi]=n.convolve(R0S[:,rgi],w,mode="same")            
    
    ED[EE > 10]=n.nan
    R0S[EE > 10]=n.nan
    EDvv[EE > 10]=n.nan    
    months=12*n.arange(365)/365
    #heights
    mm,hh=n.meshgrid(months,heights)
    plt.subplot(211)
    plt.pcolormesh(365*mm/12.0,hh,n.transpose(ED),vmin=0,vmax=30)
    plt.ylim(c.epsilon_hlimit)
    plt.xlabel("Day of year")
    plt.ylabel("Height (km)")
#    plt.title("$\\varepsilon$")    
    cb=plt.colorbar()
    cb.set_label("$\\varepsilon$ (mW/kg)")
    plt.subplot(212)
    plt.pcolormesh(365*mm/12.0,hh,n.sqrt(n.transpose(R0S)),vmin=5,vmax=25)
    plt.ylim(c.epsilon_hlimit)
    plt.xlabel("Day of year")
    plt.ylabel("Height (km)")
#    plt.title("$\sigma$")    
    cb=plt.colorbar()
    cb.set_label("$\sigma$ (m/s)")
    plt.tight_layout()
    plt.savefig("figs/epsilon.png",dpi=400)
    plt.show()


    
    if False:
        plt.subplot(222)
        plt.pcolormesh(mm,hh,n.transpose(EDvv),vmin=0,vmax=50)
        plt.ylim(c.epsilon_hlimit)
        plt.xlabel("Month")
        plt.ylabel("Height (km)")
        plt.title("$\\varepsilon_{vv}$")
        cb=plt.colorbar()
        cb.set_label("$\\varepsilon$ (mW/kg)")

        plt.subplot(223)    
        plt.pcolormesh(mm,hh,n.transpose(R0S),vmin=0,vmax=1400)
        plt.ylim(c.epsilon_hlimit)
        plt.xlabel("Month")
        plt.ylabel("Height (km)")
        cb=plt.colorbar()
        cb.set_label("$R_{LL}(0)$ (m$^2$/s$^2$)")

        plt.subplot(224)    
        plt.pcolormesh(mm,hh,n.transpose(R0Svv),vmin=0,vmax=1400)
        plt.ylim(c.epsilon_hlimit)
        plt.xlabel("Month")
        plt.ylabel("Height (km)")
        cb=plt.colorbar()
        cb.set_label("$R_{TT}(0)$ (m$^2$/s$^2$)")

        #plt.colorbar()

    #    plt.subplot(223)    

     #   plt.pcolormesh(mm,hh,n.transpose(R0S/ED),vmin=0,vmax=100,cmap="jet")
      #  plt.colorbar()
        plt.tight_layout()
        plt.show()

    # scatter plot to explore relationship between kinetic energy in
    # wind and dissipation-rate
#    plt.subplot(221)
    zero_lags=n.transpose(R0S).flatten()
    epsilon=n.transpose(ED).flatten()
    heights=hh.flatten()
    months=mm.flatten()


    if c.horizontal_fit_synoptic:
        plt.pcolormesh(n.log10(n.transpose(S)),vmin=-5,vmax=-1)
        plt.colorbar()
        plt.show()


    os.system("mkdir -p %s/%s/eps"%(c.data_directory,c.data_prefix))
    ho=h5py.File("%s/%s/eps/epsfit.h5"%(c.data_directory,c.data_prefix),"w")
    ho["R0"]=zero_lags
    ho["R_L"]=all_aacfs[:,:,0]
    ho["s_h"]=s_h
    ho["epsilon"]=epsilon
    ho["heights"]=heights
    ho["months"]=months
    ho["R0m"]=R0S
    ho["ED"]=ED    
    ho["hh"]=hh
    ho["mm"]=mm   
    ho.close()

if __name__ == "__main__":
    estimate_epsilons()
