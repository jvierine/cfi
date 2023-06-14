import h5py
import matplotlib.pyplot as plt
import numpy as n
import scipy.optimize as sio

fname="tmp.h5"
#fname="tmp_germany.h5"
h=h5py.File(fname,"r")
R0=n.sqrt(n.copy(h["R0"][()]))
eps=n.copy(h["epsilon"][()])
hgt=n.copy(h["heights"][()])
h.close()


#gidx=n.where( (R0 > n.sqrt(25)) & (R0 < n.sqrt(600)) & (eps > 1) & (eps < 50) & (hgt > 87) & (hgt < 93) )[0]
gidx=n.where( (R0 > 0) & (R0 < n.sqrt(10000)) & (eps > 0.5) & (eps < 1000) & (hgt > 0) & (hgt < 10000) )[0]

R0g=R0[gidx]
epsg=eps[gidx]

def ss(x):
    beta=x[0]
    alpha=x[1]

    return(n.sum(n.abs(beta*(R0g**alpha) - epsg)**2.0))
def ss3(x):
    beta=x[0]

    return(n.sum(n.abs(beta*(R0g**3.0) - epsg)**2.0))

def ss2(x):
    beta=x[0]

    return(n.sum(n.abs(beta*(R0g**2.0) - epsg)**2.0))

xhat=sio.fmin(ss,[1/100.0, 3.0])
xhat=sio.fmin(ss,xhat)

# linearized error estimates
# estimate the of error standard deviation
# todo: add velocity dependented error model!!!!
error_std_est=5*n.sqrt(n.mean(n.abs(xhat[0]*R0g**(xhat[1]) - epsg)**2.0))


gidx2=n.where(n.abs(xhat[0]*R0g**(xhat[1]) - epsg) < 4*error_std_est)[0]

R0g=R0g[gidx2]
epsg=epsg[gidx2]


xhat=sio.fmin(ss,[1/100.0, 3.0])
xhat=sio.fmin(ss,xhat)

#plt.plot(R0g,xhat[0]*R0g**(xhat[1]) - epsg,".")
tmps=n.linspace(0.1,30,num=100)
#plt.plot(tmps,0.5*tmps)
#plt.title(error_std_est)
#plt.show()

J = n.zeros([len(R0g),2])
xhat_dx1=n.copy(xhat)
xhat_dx2=n.copy(xhat)
xhat_dx1[0]+=1e-6
xhat_dx2[1]+=1e-6
J[:,0]=(xhat_dx1[0]*R0g**xhat_dx1[1]-xhat[0]*R0g**xhat[1])/1e-6
J[:,1]=(xhat_dx2[0]*R0g**xhat_dx2[1]-xhat[0]*R0g**xhat[1])/1e-6

J=J/error_std_est

Sp=n.linalg.inv(n.dot(n.transpose(J),J))
exp_std=2.0*n.sqrt(Sp[1,1])
print(exp_std)


xhat3=sio.fmin(ss3,[1/100.0])
xhat2=sio.fmin(ss2,[1/100.0])

plt.plot(R0[gidx],eps[gidx],".",alpha=0.1)
r0s=n.linspace(1,n.sqrt(600),num=100)
plt.plot(r0s,xhat[0]*r0s**xhat[1],label="$\\varepsilon \\propto \sigma_{u}^{%1.1f \pm %1.2f}$"%(xhat[1],exp_std))
plt.plot(r0s,xhat3[0]*r0s**3.0,label="$\\varepsilon \\propto (%1.0f)^{-1}\\sigma_{u}^{3}$"%(1/xhat[0]))
plt.plot(r0s,xhat2[0]*r0s**2.0,label="$\\varepsilon \\propto \\sigma_{u}^{2}$")
plt.legend()
plt.xlabel("$\\sigma_{u} = \\sqrt{R_{LL}(0)}$ (m/s)")
plt.ylabel("$\\varepsilon$ (mW/kg)")
plt.tight_layout()
#plt.title(xhat)
plt.show()

           
