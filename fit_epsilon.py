import numpy as n
import matplotlib.pyplot as plt
import scipy.optimize as sio





def fit_epsilon0(s,R,
                 c_k=2,
                 e0=10e-3,
                 err_std=[],
                 debug=False):
    # todo, estimate uncertainty


    def model(x):
        eps0=x[0]
        R0=x[1]
        return(R0 - 0.5*c_k*(eps0*s*1e3)**(2/3.0))

    def ss(x):
        mf=model(x)
        return(n.sum(n.abs(mf-R)**2))
        
    xhat=sio.fmin(ss,[e0,R[0]])

    resid_std_est=n.sqrt(n.mean(n.abs(model(xhat)-R)**2.0))

    # estimate errors
    xhat_dx0=n.copy(xhat)
    xhat_dx0[0]=xhat_dx0[0]+1e-6
    xhat_dx1=n.copy(xhat)
    xhat_dx1[1]=xhat_dx1[1]+1.0

    if len(err_std) == 0:
        err_std=n.zeros(len(R))
        err_std[:]=resid_std_est
        
    J=n.zeros([len(R),2])
    Sinv=n.zeros([len(R),len(R)])
    J[:,0]=(model(xhat_dx0)-model(xhat))/1e-6
    J[:,1]=(model(xhat_dx1)-model(xhat))/1.0
    for i in range(len(R)):
        Sinv[i,i]=1.0/err_std[i]**2.0
    C=n.linalg.inv(n.dot(n.dot(n.transpose(J),Sinv),J))
    eps_std=n.sqrt(C[0,0])
    r0_std=n.sqrt(C[1,1])

    if debug:
        plt.plot(s,model(xhat))
        plt.title("$\epsilon=%1.2f \pm %1.2f$ mW/kg $R_0=%1.2f \pm %1.2f$ m$^2$/s$^2$"%(xhat[0]*1e3,1e3*eps_std,xhat[1],r0_std))
        plt.plot(s,R,".")
        plt.xlabel("Horizontal separation (km)")
        plt.ylabel("Correlation function $R'_{uu}$ (m$^2$/s$^2$)")
        plt.show()

        # compensated function
        plt.plot(s,1e3*(-(2/c_k)*(R-xhat[1])*(1e3*s)**(-2/3))**(3/2),".")
        plt.plot(s,1e3*(-(2/c_k)*(R-xhat[1]))**(3/2) / (1e3*s),"x")        
        plt.axhline(xhat[0]*1e3,color="red")
        plt.title("$-\\frac{2}{C_k}[R'_{ii}(s) - R'_{ii}(0)]^{3/2} S_h^{-1}$")
        plt.xlabel("Horizontal separation (km)")
        plt.ylabel("Dissipation rate (mW/kg)")
        plt.ylim([0,100])
        plt.show()

    resid=R-model(xhat)

    
    return({"xhat":xhat,"eps_std":eps_std,"r0_std":r0_std,"resid":resid})
