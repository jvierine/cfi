import h5py
import matplotlib.pyplot as plt
import numpy as n
import scipy.optimize as sio
import scipy.interpolate as sint
from matplotlib.gridspec import GridSpec
import cfi_config as c

fname="%s/eps/epsfit.h5"%(c.data_directory)
h=h5py.File(fname,"r")
R0=n.sqrt(n.copy(h["R0"][()]))
eps=n.copy(h["epsilon"][()])
hgt=n.copy(h["heights"][()])
mnt=n.copy(h["months"][()])


plt.pcolormesh(h["mm"][()],h["hh"][()],n.transpose(n.sqrt(h["R0m"][()])**3.0/h["ED"][()]),vmin=0,vmax=1000,cmap="jet")
cb=plt.colorbar()
cb.set_label("$L_{ST} = \sigma_h^3/\varepsilon$ (km)")
plt.xlabel("Time (months)")
plt.ylabel("Height (km)")
plt.show()

h.close()
def get_std_model( r0, resid, num_step=10 ):
    r0_steps=n.linspace(n.min(r0)-1.0,n.max(r0)+1.0,num=num_step+1)
    stds=[]
    r0_model=[]
    for i in range(num_step):
        idx=n.where( (r0 > r0_steps[i]) & (r0 < r0_steps[i+1]) )[0]
        if len(idx) > 3:
            stds.append(n.std( resid[idx] ))
        else:
            stds.append(10)
        r0_model.append(0.5*(r0_steps[i]+r0_steps[i+1]))
    r0_model=n.array(r0_model)
    stds=n.array(stds)
    r0_model[0]=n.min(r0)-10
    r0_model[-1]=n.max(r0)+10
    return(sint.interp1d(r0_model, stds))
    
def fit_month( m0=0, m1=12, debug=False ):

    #gidx=n.where( (R0 > n.sqrt(25)) & (R0 < n.sqrt(600)) & (eps > 1) & (eps < 50) & (hgt > 87) & (hgt < 93) )[0]
    gidx=n.where( (R0 > 0) & (R0 < n.sqrt(10000)) & (eps > 0.5) & (eps < 1000) & (hgt > 0) & (hgt < 10000) & ((mnt > m0) & (mnt < m1)))[0]

    R0g=R0[gidx]
    epsg=eps[gidx]

    err_stds=n.zeros(len(R0g))
    err_stds[:]=1.0

    def ss(x):
        beta=x[0]
        alpha=x[1]

        return(n.sum((1/(err_stds**2.0))*n.abs(beta*(R0g**alpha) - epsg)**2.0))
    def ss3(x):
        beta=x[0]

        return(n.sum((1/(err_stds**2.0))*n.abs(beta*(R0g**3.0) - epsg)**2.0))

    def ss2(x):
        beta=x[0]

        return(n.sum((1/(err_stds**2.0))*n.abs(beta*(R0g**2.0) - epsg)**2.0))

    xhat=sio.fmin(ss,[1/100.0, 2.5])
    xhat=sio.fmin(ss,xhat)

    # linearized error estimates
    # estimate the of error standard deviation
    # todo: add velocity dependented error model!!!!

    resid = xhat[0]*R0g**xhat[1] - epsg
    stdfun=get_std_model( R0g, resid )

    err_stds=stdfun(R0g)

    #error_std_est=5*n.sqrt(n.mean(n.abs(xhat[0]*R0g**(xhat[1]) - epsg)**2.0))


    #gidx2=n.where(n.abs(xhat[0]*R0g**(xhat[1]) - epsg) < 4*error_std_est)[0]

    #R0g=R0g[gidx2]
    #epsg=epsg[gidx2]


    xhat=sio.fmin(ss,[1/100.0, 3.0])
    xhat=sio.fmin(ss,xhat)

    if debug:
        plt.plot(R0g,resid,".")
        tmps=n.linspace(0.1,30,num=100)
        plt.plot(tmps,stdfun(tmps))

        plt.show()

    J = n.zeros([len(R0g),2])
    xhat_dx1=n.copy(xhat)
    xhat_dx2=n.copy(xhat)
    xhat_dx1[0]+=1e-6
    xhat_dx2[1]+=1e-6
    J[:,0]=(xhat_dx1[0]*R0g**xhat_dx1[1]-xhat[0]*R0g**xhat[1])/1e-6
    J[:,1]=(xhat_dx2[0]*R0g**xhat_dx2[1]-xhat[0]*R0g**xhat[1])/1e-6

    for i in range(len(R0g)):
        J[i,:]=J[i,:]/err_stds[i]

    Sp=n.linalg.inv(n.dot(n.transpose(J),J))
    exp_std=2.0*n.sqrt(Sp[1,1])
    print(exp_std)


    xhat3=sio.fmin(ss3,[1/100.0])
    xhat2=sio.fmin(ss2,[1/100.0])

    return({"R0":R0[gidx],"eps":eps[gidx],"xhat":xhat,"xhat3":xhat3,"xhat2":xhat2,"exp_std":exp_std})


r=fit_month(m0=0,m1=12)
pR0=r["R0"]
peps=r["eps"]
plt.plot(pR0,peps,".",alpha=0.1)
plt.title("Full year")
r0s=n.linspace(1,n.sqrt(600),num=100)
xhat=r["xhat"]
xhat2=r["xhat2"]
xhat3=r["xhat3"]
exp_std=r["exp_std"]
plt.plot(r0s,xhat[0]*r0s**xhat[1],label="$\\varepsilon \\propto \sigma_{u}^{%1.1f \pm %1.2f}$"%(xhat[1],exp_std))
plt.plot(r0s,xhat3[0]*r0s**3.0,label="$\\varepsilon \\propto (%1.0f)^{-1}\\sigma_{u}^{3}$"%(1/xhat3[0]))
plt.plot(r0s,xhat2[0]*r0s**2.0,label="$\\varepsilon \\propto \\sigma_{u}^{2}$")
plt.legend()
plt.xlabel("$\\sigma_{u} = \\sqrt{R_{LL}(0)}$ (m/s)")
plt.ylabel("$\\varepsilon$ (mW/kg)")
plt.tight_layout()
plt.show()


fig=plt.figure(layout="constrained")
gs=GridSpec(nrows=3,
            ncols=4,
            figure=fig)

for m in range(12):
    
    row=int(n.floor(m/4))
    col=int(m-(row*4))
    ax=fig.add_subplot(gs[row,col])

    r=fit_month(m0=m,m1=m+1)
    pR0=r["R0"]
    peps=r["eps"]
    ax.plot(pR0,peps,".",alpha=0.1)
    ax.set_title(m+1)
    r0s=n.linspace(1,n.sqrt(600),num=100)
    xhat=r["xhat"]
    xhat2=r["xhat2"]
    xhat3=r["xhat3"]
    exp_std=r["exp_std"]
    ax.plot(r0s,xhat[0]*r0s**xhat[1],label="$\\varepsilon \\propto \sigma_{u}^{%1.1f \pm %1.2f}$"%(xhat[1],exp_std))
    ax.plot(r0s,xhat3[0]*r0s**3.0,label="$\\varepsilon \\propto (%1.0f)^{-1}\\sigma_{u}^{3}$"%(1/xhat3[0]))
    ax.plot(r0s,xhat2[0]*r0s**2.0,label="$\\varepsilon \\propto \\sigma_{u}^{2}$")
    ax.legend()
    ax.set_xlabel("$\\sigma_{u} = \\sqrt{R_{LL}(0)}$ (m/s)")
    ax.set_ylabel("$\\varepsilon$ (mW/kg)")
    
plt.tight_layout()
plt.show()

           
