import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv)<2:
    sys.exit("Usage: %s <file>\n"%sys.argv[0])

qn=int(sys.argv[1].split('_q')[-1].split('_m')[0])
t,terr,fsR,fsRerr,fsI,fsIerr,Navg = np.loadtxt(sys.argv[1], unpack=True)

'''
fR=fqtR[1:]
fRerr=fqtRerr[1:]

pR=fR/fqtR[0]
pRerr=pR*np.sqrt((fRerr/fR)**2 + (fqtRerr[0]/fqtR[0])**2 )
np.savetxt("phi_q%03d.dat"%qn, np.array([t[1:],pR,pRerr]).T)
'''

fig=plt.figure(figsize=(8,4))

ax1=fig.add_subplot(1,2,1)
ax1.errorbar(t[1:],fsR[1:],xerr=terr[1:],yerr=fsRerr[1:],fmt="b.-")
ax1.set_xscale('log')
ax1.tick_params(which='both', direction='in')
ax1.set_xlabel("t / ps")
ax1.set_ylabel("Re[F$_s$(q,t)]")

ax3=fig.add_subplot(1,2,2)
ax3.axhline(0,c='k',ls='--')
ax3.errorbar(t[1:],fsI[1:],xerr=terr[1:],yerr=fsIerr[1:],fmt="k.-")
ax3.set_xscale('log')
ax3.tick_params(which='both', direction='in')
ax3.set_xlabel("t / ps")
ax3.set_ylabel("Im[F$_s$(q,t)]")

plt.title("Average over %d-%d data points"%(min(Navg),max(Navg)))
plt.tight_layout()
#fig.savefig("fself_ave_q%03d.eps"%qn)

plt.show()
