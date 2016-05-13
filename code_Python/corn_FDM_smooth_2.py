# -*- coding: utf-8 -*-
"""
digital option(cash-or-nothing) version 2 smoothing
@author: Minhyun Yoo
"""
from time import time
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, exp
from scipy import stats

# thomas algorithm for tridiagonal matrix
def thomas(alpha, beta, gamma, f):
    n = len(f)
    x = np.zeros(n)
    
    # To treat the shallow copy, map is used.
    # refence : http://fab.cba.mit.edu/classes/864.14/students/Langford_Will/5_fd_pdes/
    ac, bc, gc, fc = map(np.array, (alpha, beta, gamma, f)) # copy the array
    
    for i in xrange(1,n):
        mult = ac[i] / bc[i-1]
        bc[i] = bc[i] - mult*gc[i-1]
        fc[i] = fc[i] - mult*fc[i-1]
    x[n-1] = fc[n-1] / bc[n-1]
    
    for i in xrange(n-2,-1,-1):
        x[i] = (fc[i] - gc[i] * x[i+1]) / bc[i]
    
    del ac, bc, gc, fc
    return x

timer = time()
Nt = 100; # # of time steps
T = 1. / 365.; # Near the maturity
dt = T / Nt; # For around the maturity
E = 100.; # strike price
L = 300.; # sufficiently large value(computational domain)
sig = 0.5; # volatility
r = 0.03; # riskless interest rate
dx = 0.1; # spatial step

# grid construct - nonuniform
x = np.concatenate([np.array([0]), np.arange(0.5*dx,L-0.5*dx,dx)]);
Nx = len(x);
Nt = 100;
h = np.diff(x);
h = np.concatenate([h, np.array([h[-1]])]);

# digital option(cach-or-nothing)
cash = 100;
u = np.zeros((Nx, Nt+1));
###### version 2 ######
equd = 0.5; # can be modified
rev = (x - E) / equd;

y2 = np.zeros(Nx);
for i in xrange(Nx):
    if ((rev[i] >= -1.0) & (rev[i] <= 0.0)):
        y2[i] = 0.5*(1.0 + rev[i])**2;
    elif ((rev[i] >= 0.0) & (rev[i] <= 1.0)):
        y2[i] = 1.0 - 0.5*(1.0 - rev[i])**2;
    elif (rev[i] < -1.0):
        y2[i] = 0.0;  
    else:
        y2[i] = 1.0;
u[:, 0] = cash * y2;
########################

# tridiagonal matrix
f = np.zeros(Nx-1); d = np.zeros(Nx-1); 
a = np.zeros(Nx-1); c = np.zeros(Nx-1);
for i in range(1, Nx):
    d[i-1] = 1.0 + dt * (((sig*x[i])**2 - r*x[i]*(h[i]-h[i-1])) / (h[i-1]*h[i]) + r );
    c[i-1] = dt * ( -(sig*x[i])**2 - r*x[i]*h[i-1]) / (h[i]*(h[i]+h[i-1]));
    a[i-1] = dt * ( -(sig*x[i])**2 + r*x[i]*h[i]) / (h[i-1]*(h[i]+h[i-1]));
# linear boundary condition
d[Nx-2] = d[Nx-2] + 2.0*c[Nx-2];
a[Nx-2] = a[Nx-2] - c[Nx-2];

# time loop
for n in range(Nt):
    b = u[1:Nx, n];
    u[1:Nx, n+1] = thomas(a,d,c,b);
    
timer = time() - timer;
print "Duration in Seconds %7.5f" % timer

# exact option price
d1 =  (np.log(x/E) + (r + 0.5*sig**2)*T) / (sig*sqrt(T));
d2 = d1 - (sig*sqrt(T));
exc = cash * exp(-r * T) * stats.norm.cdf(d2);

# error
bc = 0.8; ec = 1.2;
# find index 
# reference : http://stackoverflow.com/a/25032853/5861666
bidx =  np.searchsorted(x, bc * E);
eidx = np.searchsorted(x, ec * E)-1;

RMSE = np.sqrt( np.mean( (u[bidx:eidx,-1] - exc[bidx:eidx])**2 ));
maxerr = np.max(np.abs(u[:,-1] - exc));

print "RMSE : %.8f" % RMSE;
print "maxerr : %.8f" % maxerr;

# plot
# secondary axis
# reference : http://stackoverflow.com/a/5487005/5861666
fig = plt.figure(figsize = (16,10), dpi = 100)
plt.rcParams.update({'font.size': 18})
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
pl1 = ax1.plot(x, u[:,0], c='k', label='Payoff')
pl2 = ax1.plot(x, u[:,Nt], c='b', marker='.', label='FDM solution')
pl3 = ax1.plot(x, exc, c='r', label='Exact')
pl4 = ax2.plot(x, np.abs(u[:,-1] - exc), c='r', linestyle ='--', label='Error')
ax1.grid(True)
ax1.set_xlim(90, 110)
ax1.set_xlabel('$x$')
ax1.set_ylabel('$u(x,t)$')
ax2.set_ylabel('$error$')
plt.title('corn_FDM_smooth_2.py')

lns = pl1+pl2+pl3+pl4;
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=0)
