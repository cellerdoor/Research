import numpy as np
import matplotlib.pyplot as plt
import random
from math import *
from input import *

xs = np.linspace(-L/2., L/2., L/dx)
Nu = len(xs)
A0 = kon/gamma
dA = 0.1*A0
Am = A0 + dA*2.0*(0.5-np.random.randn(Nu))
t = 0.0
ub = 1.0
ubs = [ub]
dt = 0.05
ub2 = ub - .5
ubs2 = [ub2]

ts = [t]
Ams = [list(Am)]
am = list(Am)
for i in range(10000):

    #force calculations
    fx = pi*np.trapz([exp(-0.5*(xsi-ub)**2.0)*(xsi-ub)*Ami for xsi,Ami in zip(xs,Am)], xs )
    fx2 = pi*np.trapz([exp(-0.5*(xsi-ub2)**2.0)*(xsi-ub2)*Ami for xsi,Ami in zip(xs,Am)], xs )
    v0s = [exp(-0.5*(xsi-ub)**2.0/c/c) for xsi in xs]
    v0s2 = [exp(-0.5*(xsi-ub2)**2.0/c/c) for xsi in xs]
    Am = [Ami+dt*(kon-Ami-gamma*v0i*Ami) for Ami,v0i in zip(Am,v0s)]
    Am = [Ami+dt*(kon-Ami-gamma*v0i*Ami) for Ami,v0i in zip(Am,v0s2)]
    
    # update position
    ub += fx*dt
    ubs.append(ub)
    ub2 += fx2*dt
    ubs2.append(ub2)
    t += dt
    ts.append(t)


plt.plot(xs,Am)
plt.show()
'''
ub3 = ub2+.1
ubs3 = [ub3]
ub4 = ub - .1
ubs4 = [ub4]
for i in range(10000):
    #force calculations
    fx = pi*np.trapz( exp(-0.5*(xs-ub)**2.0)*(xs-ub)*Am, xs )
    fx2 = pi*np.trapz( exp(-0.5*(xs-ub2)**2.0)*(xs-ub2)*Am, xs )
    fx3 = pi*np.trapz( exp(-0.5*(xs-ub3)**2.0)*(xs-ub3)*Am, xs )
    fx4 = pi*np.trapz( exp(-0.5*(xs-ub4)**2.0)*(xs-ub4)*Am, xs )

    v0s = exp(-0.5*(xs-ub)**2.0/c/c)
    v0s2 = exp(-0.5*(xs-ub2)**2.0/c/c)
    v0s3 = exp(-0.5*(xs-ub3)**2.0/c/c)
    v0s4 = exp(-0.5*(xs-ub4)**2.0/c/c)

    Am += kon*dt - gamma*v0s*Am*dt - gamma*v0s2*Am*dt - gamma*v0s3*Am*dt - gamma*v0s4*Am*dt - Am*dt

    # update position
    ub += fx*dt 
    ubs.append(ub) 
    ub2 += fx2*dt 
    ubs2.append(ub2)
    ub3 += fx3*dt 
    ubs3.append(ub3)
    ub4 += fx4*dt 
    ubs4.append(ub4)

    t += dt
    ts.append(t)
'''
plt.plot(ts,ubs)
plt.plot(ts,ubs2)
#plt.plot(ts,ubs2)
#plt.plot(ts[-len(ubs3):],ubs3)
#plt.plot(ts[-len(ubs4):],ubs4)
plt.show()


