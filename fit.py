import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import gamma as Gamma
import sys
import math

gamma=2 # exponential or gaussian decay? 1 : 2

def strexp(t,tau=1,beta=1,a=.5):
    return a*np.exp(-(t/tau)**beta)

def phonons(t,fq=.5,tau=1,omega=1):
    return strexp(t,tau,gamma,1-fq)*np.cos(omega*t) + fq

def fqt(t,tauA=1,beta=1, fq=.5, tau=1,omega=1):
    return strexp(t,tau,gamma,1-fq)*np.cos(omega*t) + strexp(t,tauA,beta,fq)

fitmode=-1
fitfunc=fqt
fitstr="phonons+decay fit"
params=["tauA","beta","fq","tau","omega"]
p0=    [   1e2,   0.6, 0.5,    1,    3]
lowb=np.array(
       [   0.1, 1e-2,    0, 1e-2,  1e-3]
)
uppb=np.array(
       [np.inf,  1e2,    1,  1e3,   1e2]
)

if len(sys.argv)<2 or len(sys.argv)>6:
    sys.exit("USAGE: %s <phi-datafile> [gamma] [t1] [t2] [fitting-mode]\n\n Fits normalized ISF(q,t) at fixed wavevector q with a short-time partially decaying oscillation (phonons) and a long-time stretched exponential relaxation).\n\nARGUMENTS:\n -phi-datafile: path to a file containing t,phi,phi-error in columns.\n -gamma: exponential or gaussian decay? 1 : 2 (default %d)\n -t1,t2 (optional): fitted time interval, in picoseconds (default: total data)\n -Fitting mode (optional): 0: only phonons; 1: only relaxation; (default: sum of them).\n"%(sys.argv[0],gamma))

path,tail=sys.argv[1].split("phi_q")
qn=int(tail.split("_m")[0])#int(sys.argv[1].split('.')[0].split('q')[-1])
M=int( tail.split("_m")[-1].split("_nat")[0] )
NAT=int( tail.split("_nat")[-1].split(".")[0] )
t,phi,phierr = np.loadtxt(sys.argv[1], unpack=True)
t=t[1:]
phi=phi[1:]
phierr=phierr[1:]
t1=0
t2=t[-1]
slice=np.ones(len(t), dtype=bool)
if len(sys.argv)>2:
    gamma=int(sys.argv[2])
    if len(sys.argv)>3:
        t1=float(sys.argv[3])
        slice = (t>=t1)
        if len(sys.argv)>4:
            t2=float(sys.argv[4])
            slice = (t>=t1)*(t<=t2)
            if len(sys.argv)>5:
                fitmode=int(sys.argv[5])
outfile = f"{path}fit_q{qn:03d}_m{M:03d}_nat{NAT:d}_g{gamma:d}"

if fitmode==0:
    fitfunc=phonons
    fitstr="phonons fit"
    params=params[2:]
    p0=p0[2:]
    lowb=lowb[2:]
    uppb=uppb[2:]
elif fitmode==1:
    fitfunc=strexp
    fitstr="decay fit"
    params=params[:3]
    p0=p0[:3]
    lowb=lowb[:3]
    uppb=uppb[:3]

bounds=(lowb,uppb)
x=t[slice]
y=phi[slice]
yerr=phierr[slice]


xl=np.logspace(np.log10(t[0]),np.log10(t[-1]),num=1000,base=10.0)
try:
    popt,pcov = curve_fit(fitfunc, x,y,sigma=yerr, p0=p0,bounds=bounds, maxfev=10000)
except RuntimeError:
    print("RuntimeError: curve_fit of phi(t;q=%d) failed."%qn)
    sys.exit(1)
f=open(outfile+".log","w")
print("#Time interval: %f to %f ps"%(x[0],x[-1]))
print("")
f.write("#Time interval: %f to %f ps\n"%(x[0],x[-1]))
for i in range(len(popt)):
    print("%10s=%6.8f +- %6.8f"%(params[i],popt[i],np.sqrt(pcov[i][i])))
    f.write("%10s=%6.8f +- %6.8f\n"%(params[i],popt[i],np.sqrt(pcov[i][i])))
print("")
### now calculate the average t along a stretched exponential distribution
### it's a better estimate for tau-alpha (?)
samples=np.random.multivariate_normal(popt[:2],pcov[:2,:2],size=1000)
tmedio=popt[0]*Gamma(1/popt[1])/popt[1]
tmedioerr=(samples[:,0]*Gamma(1/samples[:,1])/samples[:,1]).std() #estimate the error as the std.dev. of 1000 values extracted from multivariate gaussian
print("%10s=%6.8f +- %6.8f"%("tmedio",tmedio,tmedioerr))
f.write("%10s=%6.8f +- %6.8f\n"%("tmedio",tmedio,tmedioerr))
Y = fitfunc(x,*popt)
chi2 = np.sum( ((y-Y)/yerr)**2 )/ (len(y)-1)
print("%10s=%6.8f"%("chi2",chi2))
f.write("%10s=%6.8f\n"%("chi2",chi2))
###
f.close()

yl=fitfunc(xl,*popt)
np.savetxt(outfile+".dat", np.array([xl,yl,np.zeros(len(xl))]).T)

'''# plot initial guess 
fig0=plt.figure()
ax=fig0.add_subplot(1,1,1)
ax.plot(xl,yl)
ax.set_xscale('log')
ax.set_xlabel('t / ps')
ax.set_ylabel('$\phi$(q,t)')
ax.set_title('Initial guess')
'''

# plot data and optimized fitting function
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
#ax.axhline(0,c='k',ls='--')
ax.errorbar(t,phi,phierr,fmt="k.", label="data")
ax.errorbar(x,y,yerr,fmt="b.", label="fitted interval")#[%.1e, %.1e] ps"%(t1,t2))
ax.plot(xl,yl,"r",label=fitstr)
ax.set_xscale('log')
ax.tick_params(which='both', direction='in')
ax.set_xlabel("t / ps")
ax.set_ylabel("$\phi(q,t)$")
ax.legend()

#fig.savefig("%sfit_q%03d.eps"%(path,qn))
plt.show()
